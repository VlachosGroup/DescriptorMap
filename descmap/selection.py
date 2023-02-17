import pandas as pd
import numpy as np
from itertools import combinations
from sklearn.decomposition import PCA
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.model_selection import LeaveOneOut
from sklearn.linear_model import LinearRegression

def _get_energies(model, T):
    """Helper method to calculate adsorption energies of species using PCA model

    Parameters
    ----------
        model : sklearn.decomposition.PCA object
            Trained PCA model
        T : (n_samples, n_components)
            Principal components matrix
    """
    return model.inverse_transform(T)

def _get_robustness(model, X, energies):
    """Helper function to calculate robustness metric from species energies using PCA model
    
    Parameters
    ----------
        model : sklearn.decomposition.PCA.object
            Trained PCA model
        X : pandas DataFrame object
            (m_surfaces, n_adsorbates). Adsorption energies of species read from excel
        energies: Array of float
            Calculated adsorption energies of species using PCA model
    Returns
    -------
        robustness_vec : pandas Series object
    """
    rmse_vec = np.zeros(len(X.columns))
    robustness_vec = np.zeros(len(X.columns))
    for i, ads in enumerate(X.columns):
        rmse_vec[i] = mean_squared_error(X.iloc[:, i],
                                         energies[:, i],
                                         squared=False)
        robustness = 0.
        for weight, eigenvec in zip(model.components_[:, i], model.explained_variance_):
            robustness += eigenvec*(weight)**2
        robustness_vec[i] = np.sqrt(robustness)/rmse_vec[i]
    return pd.Series(robustness_vec, index=X.columns)

def get_component_number(X, var_explained):
    """Calculate the number of PCA components to reach the explained variance
    
    Parameters
    ----------
        X : pandas DataFrame object
            (m_surfaces, n_adsorbates). Adsorption energies of species read from excel
        var_explained: float
            Defined thresold for explained variance
    
    Returns
    -------
        n_components : int
            Required number of components
    """   
    # Set up PCA inputs
    model = PCA()
    Y = pd.DataFrame(model.fit_transform(X), index=X.index)
    
    # Calculate the number of components required
    norm_cum_variance = np.cumsum(model.explained_variance_/np.sum(model.explained_variance_))
    norm_cum_variance_cutoff = norm_cum_variance[norm_cum_variance < var_explained]
    n_components = len(norm_cum_variance_cutoff) + 1
    return n_components

def get_combination_score(X, n_components):
    """Calculate the metrics for each combination of descriptor set
    
    Parameters
    ----------
        X : pandas DataFrame object
            (m_surfaces, n_adsorbates). Adsorption energies of species read from excel
        n_components: int
            Number of principal components for PCA model
    
    Returns
    -------
        comb_df : pandas Series object
            Fitness score for each descriptor combination
        complete_df : pandas DataFrame object
            Calculated metrics for each descriptor combination: inv_cond_num, robustness_vec, fitness    
    """   
    # Calculate adsorption energies using PCA model
    model_var_cutoff = PCA(n_components=n_components)
    Y_var_cutoff = model_var_cutoff.fit_transform(X)
    components_var_cutoff = pd.DataFrame(model_var_cutoff.components_, columns=X.columns)
    energies = _get_energies(model_var_cutoff, Y_var_cutoff)
    
    # Calculate robustness vector
    robustness_vec = _get_robustness(model_var_cutoff, X, energies)
    norm_robustness_vec = (robustness_vec - min(robustness_vec))/(max(robustness_vec) - min(robustness_vec))
    
    # Generate descriptor combination dataframe
    species_comb_list = []
    fitness_comb_list = []
    inv_cond_num_list = []
    robustness_vec_comb_list = []
        
    for k, species_comb in enumerate(combinations(X.columns, r=n_components)):
        basis_func = np.array([components_var_cutoff.loc[:, species] for species in species_comb])
        inv_cond_num = 1./np.linalg.cond(basis_func)
        
        robust_weight = 0.5/n_components
        robustness_vec_comb = np.array([norm_robustness_vec.loc[species] for species in species_comb])
        
        fitness = 0.5*inv_cond_num + np.sum(robust_weight*robustness_vec_comb)
        fitness_comb_list.append(fitness)
        species_comb_list.append(', '.join(species_comb))
        inv_cond_num_list.append(inv_cond_num)
        robustness_vec_comb_list.append(robustness_vec_comb)
    
    comb_data = {'inv_cond_num': inv_cond_num_list,
                 'robustness_vec': robustness_vec_comb_list,
                 'fitness': fitness_comb_list}
    complete_df = pd.DataFrame(comb_data, index=species_comb_list)
    comb_df = pd.Series(fitness_comb_list, index=species_comb_list)
    return comb_df, complete_df

def get_best_combination(X, n_max_combs, comb_df):
    """Get the best descriptor combination with the smallest mean squared error (MSE)
    
    Parameters
    ----------
        X : pandas DataFrame object
            (m_surfaces, n_adsorbates). Adsorption energies of species read from excel
        n_max_combs: int
            Maximum number of descriptor combinations to consider
        comb_df: pandas DataFrame object
            Fitness score for each descriptor combination
    
    Returns
    -------
        comb_err : pandas DataFrame object
            MSE and MAE for each descriptor combination
        best_comb : pandas DataFrame object
            Best descriptor combination with corresponding MSE and MAE    
    """
    comb_df = comb_df.nlargest(n=n_max_combs, keep='all').append(pd.Series([comb_df.loc['H, C']], index=['H, C']))
    comb_mses = np.zeros(n_max_combs+1)
    comb_maes = np.zeros(n_max_combs+1)
    for i, species_comb_str in enumerate(comb_df.index):
        species_comb = tuple(species_comb_str.split(', '))
        y_species = list(set(X.columns)-set(species_comb))
        X_pred = X.loc[:, species_comb]
        Y = X.loc[:, y_species]
        cv_splitter = LeaveOneOut()
        lin_reg_mses = []
        lin_reg_maes = []
        for train_i, val_i in cv_splitter.split(X_pred):
            X_train, X_val = X_pred.iloc[train_i, :], X_pred.iloc[val_i, :]
            Y_train, Y_val = Y.iloc[train_i, :], Y.iloc[val_i, :]

            lin_reg_model = LinearRegression(fit_intercept=True)
            lin_reg_model.fit(X=X_train, y=Y_train)

            Y_val_model = lin_reg_model.predict(X_val)
            lin_reg_mse = mean_squared_error(Y_val_model, Y_val)
            lin_reg_mses.append(lin_reg_mse)

            lin_reg_mae = mean_absolute_error(Y_val_model, Y_val)
            lin_reg_maes.append(lin_reg_mae)
        comb_mses[i] = np.mean(lin_reg_mses)
        comb_maes[i] = np.mean(lin_reg_maes)
    comb_err = pd.DataFrame(np.array([comb_mses, comb_maes]).T,
                            index=comb_df.index,
                            columns=['MSE', 'MAE'])
    best_comb = comb_err.nsmallest(n=1, columns='MSE')
    return comb_err, best_comb

def get_lsr_results(X, best_comb):
    """Get the LSR slopes and intercepts using the best descriptor combination
    
    Parameters
    ----------
        X : pandas DataFrame object
            (m_surfaces, n_adsorbates). Adsorption energies of species read from excel
        best_comb : pandas DataFrame object
            Best descriptor combination with smallest MSE
    
    Returns
    -------
        final_df : pandas DataFrame object
            Intercepts and slopes for each species
    """
    final_lin_reg_model = LinearRegression(fit_intercept=True)
    best_comb_tuple = tuple(best_comb.index[0].split(', '))
    y_species = list(set(X.columns)-set(best_comb_tuple))
    X_final = X.loc[:, best_comb_tuple]
    Y_final = X.loc[:, y_species]

    final_lin_reg_model.fit(X=X_final, y=Y_final)
    Y_final_pred = final_lin_reg_model.predict(X_final)
    final_data = {'{} intercept'.format(best_comb_tuple[0]): final_lin_reg_model.coef_[:, 0],
                  '{} intercept'.format(best_comb_tuple[1]): final_lin_reg_model.coef_[:, 1],
                  'intercept': final_lin_reg_model.intercept_}
    final_df = pd.DataFrame(data = final_data, index = y_species)
    return final_df

def write_lsr_results(final_df):
    """Write intercepts and slopes to excel for future use
    
    Parameters
    ----------
        final_df : pandas DataFrame object
            Intercepts and slopes for each species
    """
    final_df.to_excel('lsr_results.xlsx')