import numpy as np
import xgboost as xgb

xgb_oracle = None
    
def xgb_oracle_predict(triplet, cell_df, m_i, n_i, fp_i, fn_i):
    """Uses XGBoost-based oracle to predict split-away cell

    Args:
        triplet ([string]): cell triplet query
        cell_df (pandas.DataFrame): genotype matrix
        m_i (int): mutation count
        n_i (int): cell count
        fp_i (float): false positive rate
        fn_i (float): false negative rate

    Returns:
        string: split-away cell
    """
    global xgb_oracle
    # cache oracle to avoid repeatedly loading
    if not xgb_oracle:
        
        xgb_oracle_fp = f"sample-oracles/xgboost-oracle-simNo_{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}}-m_{m_i}-n_{n_i}-fp_{fp_i}-fn_{fn_i}.ubj"
        xgb_oracle = xgb.Booster()
        xgb_oracle.load_model(xgb_oracle_fp)
    # prepare input for oracle format
    cat_geno  = np.concat(cell_df.loc[triplet].values).reshape((1, -1))
    test_inpt = xgb.DMatrix(cat_geno)
    idx = int(xgb_oracle.predict(test_inpt)[0])
    return triplet[idx]
    
    
    