import numpy as np

def get_lambda_nu():
    
    """
    Return a vector of precomputed lambda_nu values. These values correspond to np = 1:13.
    
    Output:
        lambda_vec (float array) - array containing values of lambda_nu. The nth element contains lambda_nu
                                   evaluated for np = n+1
    
    """
    
    lamb1 = 14.6666666; 
    lamb2 = 45.0520048;
    lamb3 = 106.3140869;
    lamb4 = 229.0403514;
    lamb5 = 474.5838228;
    lamb6 = 965.7138504;
    lamb7 = 1947.9948706;
    lamb8 = 3912.5672515;
    lamb9 = 7841.71714841;
    lamb10 = 15700;
    lamb11 = 31416.6;
    lamb12 = 62849.8;
    lamb13 = 125716;
    lamb14 = 251432;
    lamb15 = 502864;
    lamb16 = 1005728;
    lamb17 = 2011456;
    lamb18 = 4022912;
    lamb19 = 8045824;
    lamb20 = 16091648;
    lamb21 = 32183296;
    lamb22 = 64366592;
    
    lambda_vec = [lamb1,lamb2,lamb3,lamb4,lamb5,lamb6,lamb7,lamb8,lamb9,lamb10,lamb11,lamb12,lamb13,lamb14,lamb15,lamb16,lamb17,lamb18,lamb19,lamb20,lamb21,lamb22]
    
    return lambda_vec
