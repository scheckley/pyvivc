# Python port of Rivivc numerical convolution by Aleksander Mendyk, Sebastian Polak 
# https://cran.r-project.org/web/packages/Rivivc/index.html
# Version 1.1
# Stephen Checkley, September 2021

from scipy.interpolate import pchip_interpolate as pchip
from scipy import stats
from scipy.optimize import minimize, dual_annealing
import pandas as pd
import numpy as np
from newrange import *

###############################
# Numerical convolution       #
###############################

def NumConv(impulse, input, conv_timescale=None, explicit_interpolation=1000):

    accuracy = explicit_interpolation

    input_orig = input

    impulse_orig = impulse

    # check convolution timescale
    if conv_timescale is None:
        time_orig = input_orig['time']

    else:
        if len(conv_timescale) > 2:
            time_orig = conv_timescale
        else:
            time_orig = crange(conv_timescale.iloc[0], conv_timescale.iloc[1], (conv_timescale.iloc[1]-conv_timescale.iloc[0])/accuracy)

    time_imp_orig = impulse_orig['time']
    lk_row_orig = len(time_orig)-1

    time = crange(time_orig.iloc[0], time_orig.iloc[lk_row_orig], (time_orig.iloc[lk_row_orig]-time_orig.iloc[0])/accuracy)

    lk_row1 = len(time)-1

    input1 = pchip(time_orig, input_orig['C'], time)
    impulse1 = pchip(time_imp_orig, impulse_orig['C'], time)

    def PKconvolution(input, impulse, lk_row):
        convolution = np.zeros(lk_row)
                        
        for i in range(0,lk_row-1):
            pom = 0
            for k in range(0,i):
                input1 = input[i-k+1]
                input2 = 0

                if i != k:
                    input2 = input[i-k]

                if k==1:
                    pom = pom+(impulse[k])/2*(input1-input2)
                else:
                    pom = pom+(impulse[k]+impulse[k-1])/2*(input1-input2)
            
                convolution[i] = pom

                                
        return convolution

    convol_res = PKconvolution(input1, impulse1, lk_row1)

    out_values_time = np.zeros(len(convol_res)-1)
    out_values_res = np.zeros(len(convol_res)-1)

    for i in range(len(convol_res)-1):
        out_values_time[i] = time[i]
        out_values_res[i] =  convol_res[i]

    out = pd.DataFrame()
    out['time'] = out_values_time
    out['convolution'] = out_values_res

    return out

###############################
# Numerical deconvolution     #
###############################

def NumDeconv(impulse, response, dose_iv=None, dose_po=None, deconv_timescale=None, explicit_interpolation=15, implicit_interpolation=10, maxit_optim=200, global_optim=0):

    iv_multiplication_factor = 1

    if dose_iv and dose_po is not None:
        iv_multiplication_factor = dose_po/dose_iv

    accuracy = explicit_interpolation
    multipl_2 = implicit_interpolation
    maxit_optim = maxit_optim

    # explicit interpolation
    impulse_orig = impulse
    resp_orig = response
    time_resp_orig = resp_orig['time']
    time_imp_orig = impulse_orig['time']

    if iv_multiplication_factor != 1:
        for i in range(len(impulse_orig)):
            impulse_orig['C'][i] = impulse_orig['C'][i] * iv_multiplication_factor

    # check deconvolution timescale
    if deconv_timescale is None:
        time_orig = time_imp_orig
    else:
        if len(deconv_timescale) > 2:
            time_orig = deconv_timescale
        else:
            time_orig = crange(deconv_timescale.iloc[0], deconv_timescale.iloc[1], (deconv_timescale.iloc[1]-deconv_timescale.iloc[0])/accuracy)

    lk_row_orig = len(time_orig)-1

    time = crange(time_orig.iloc[0], time_orig.iloc[lk_row_orig], (time_orig.iloc[lk_row_orig]-time_orig.iloc[0])/accuracy)
    
    lk_row1 = len(time)-1
    
    impulse_interp = pchip(time_imp_orig, impulse_orig['C'], time)
    resp_interp = pchip(time_resp_orig, resp_orig['C'], time)

    time_2 = crange(time_orig.iloc[0], time_orig.iloc[lk_row_orig], ((time_orig.iloc[lk_row_orig] - time_orig.iloc[0])/(multipl_2*accuracy)))

    resp_interp_2 = pchip(time_resp_orig, resp_orig['C'], time_2)
    impulse_interp_2 = pchip(time_imp_orig, impulse_orig['C'], time_2)
    lk_row_2 = len(time_2)
    
    #################################################################

    def error_function(input_data):
        def MSE(vect1, vect2):
            result = 0
            rows_no = len(vect1)

            for i in range(0,rows_no-1):
                result = result + np.square(np.subtract(vect1[i], vect2[i])).mean() # MSE
                #result = result + np.sum((vect1[i] - vect2[i])**2) # LSE

            return result
        
        def PKconvolution(input, impulse, lk_row):
            convolution = np.zeros(lk_row)
                        
            for i in range(0,lk_row-1):
                pom = 0
                for k in range(0,i):
                    input1 = input[i-k+1]
                    input2 = 0

                    if i != k:
                        input2 = input[i-k]

                    if k==1:
                        pom = pom+(impulse[k])/2*(input1-input2)
                    else:
                        pom = pom+(impulse[k]+impulse[k-1])/2*(input1-input2)
            
                    convolution[i] = pom

                                
            return convolution

        error = 1000

        input_data_2 = pchip(time, input_data, time_2)

        convol_res = PKconvolution(input_data_2, impulse_interp_2, lk_row_2)

        error = MSE(convol_res, resp_interp_2)
        
        return error

    # initialize the input vector values


    input_discovered = np.ones(lk_row1+1)
    
    #multipl = 1/lk_row1

    #for i in list(range(0,lk_row1+1)):
    #    input_discovered[i] = i*multipl



    # optimize
    
    print('running optimization...')

    bounds = list(zip(np.zeros(len(input_discovered)), np.ones(len(input_discovered)*100)))
    if global_optim==1:
        fit = dual_annealing(func=error_function, bounds=bounds, seed=6174)
    else:
        fit = minimize(fun=error_function, x0=input_discovered, bounds=bounds, method="L-BFGS-B", options={'maxiter':maxit_optim})
    

    input_discovered = fit.x

    interpolated = pchip(time, input_discovered, time_orig)
    interpolated_2 = pchip(time, input_discovered, time_2)
    
    out_values_time = np.zeros(lk_row_orig)
    out_values_interpolated = np.zeros(lk_row_orig)
    out_values_1_time = np.zeros(lk_row1)
    out_values_1_interpolated = np.zeros(lk_row1)
    out_values_2_time = np.zeros(lk_row_2)
    out_values_2_interpolated = np.zeros(lk_row_2)

    F_time = np.zeros(lk_row_2)
    F_out = np.zeros(lk_row_2)

    # orig_scale

    for i in range(lk_row_orig):
        out_values_time[i] = time_orig.iat[i]
        out_values_interpolated[i] =  interpolated[i]

    # epxlicit interpolation
    for i in range(lk_row1):
        out_values_1_time[i] = time[i]
        out_values_1_interpolated[i] = input_discovered[i]

    # implicit interpolation
    for i in range(lk_row_2):
        out_values_2_time[i] = time_2[i]
        out_values_2_interpolated[i] = interpolated_2[i]


    par_out = pd.DataFrame()
    par_out['time'] = out_values_time
    par_out['par'] = out_values_interpolated
    
    explicit_out = pd.DataFrame()
    explicit_out['time'] = out_values_1_time
    explicit_out['par'] = out_values_1_interpolated

    inv_explicit_out = pd.DataFrame()
    inv_explicit_out['time'] = explicit_out['time']
    inv_explicit_out['par'] = 1/explicit_out['par']

    implicit_out = pd.DataFrame()
    implicit_out['time'] = out_values_2_time
    implicit_out['par'] = out_values_2_interpolated

    Ft = np.gradient(implicit_out['par'], implicit_out['time']) # Calculate the gradient of the deconvolution, the slope of which is equivalent to Ka
    
    Ft_dat = pd.DataFrame()
    Ft_dat['time'] = implicit_out['time']
    Ft_dat['Ft'] = Ft

    return [par_out, explicit_out, implicit_out, Ft_dat]


###############################
# IVIVC A                     #
###############################                 

def pyivivc(known_dat, impulse, second_profile, dose_iv=None, dose_po=None, explicit_interpolation=20, implicit_interpolation=10, maxit_optim=200):
    x = known_dat['C']
    known_time = known_dat['time']

    wynik = NumDeconv(impulse, second_profile, dose_iv=dose_iv, dose_po=dose_po, deconv_timescale=known_time, explicit_interpolation=explicit_interpolation, implicit_interpolation=implicit_interpolation, maxit_optim=maxit_optim)
    y = wynik[0]['par']

    regr = stats.linregress(x[0:len(x)-1],y)

    out_regression = regr
    out_numeric = wynik[0]

    return [out_regression, out_numeric]
