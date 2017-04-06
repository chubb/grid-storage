import numpy as np

def PV(rate, payments):
    """Find the present value for each payment in an array of periodic payments.
    
    Args:
        rate: The discount rate for the period.
        payments: A numpy array of periodic payments.

    Returns:
        A numpy array of the present value of each payment.
    
    Example:  
        rate = 0.04/12.0  # 4% annual interest rate compounded monthly.
        payments = np.ones(12*2)*100.0  # Two years of $100 payments.
        pv_payments = PV(rate, payments)
    """
    pv_payments = np.zeros(np.size(payments))
    discount = 1.0
    for i in range(len(payments)):
        discount *= 1.0 + rate
        pv_payments[i] = payments[i]/discount
    return pv_payments