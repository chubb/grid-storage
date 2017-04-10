import numpy as np


def UpdateBatteryCapacity(capacity_0, capacity, energy_stored):
    """Find the battery capacity as a function of previous use.
    
    Args:
       capacity_0: The initial battery capacity.
       capacity: The current battery capacity.
       energy_stored: Energy stored during this cycle.
       
    Returns:
       The updated battery capacity.
    """
    assert(energy_stored <= capacity_0)

    dod = energy_stored/capacity_0

    (num_cycles, capacity_at_eol) = BatteryLife60(dod)

    # Find the percentage reduction in capacity.
    reduction_in_capacity = (1.0 - capacity_at_eol) * 1.0/num_cycles

    return capacity - reduction_in_capacity * capacity_0


def BatteryLife60(dod):
    """ Find the battery life in number of cycles given a depth of discharge.

    Args:
        dod: Depth of Discharge
    
    Returns:
        num_cycles: The number of cycles to reach 60% of original capacity.
    
    References:
        li_ion_battery_life__TechnicalSheet_en_0514_Protected.pdf
    """
    assert(np.all(0.0 <= dod) and np.all(dod <= 1.0))

    # This data is for EoL defined as 60%
    capacity_at_eol = 0.6

    dods = np.array([0.01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0])
    num_cycles = np.array([1.00E+07, 1.00E+06, 3.00E+05, 1.00E+05, 4.00E+04,
                           2.10E+04, 1.30E+04, 1.00E+04, 8.00E+03, 7.00E+03,
                           6.50E+03, 6.00E+03])
    
    # Most accounts from those who have worked in Li Ion say that
    # actual cycles are much less.
    sales_pitch_factor = 0.5

    cycles_at_eol = LogLogInterp(dod, dods, num_cycles*sales_pitch_factor)

    return (cycles_at_eol, capacity_at_eol)


def LogLogInterp(x, xs, ys):
    """Interpolation well suited for data that is linear on a log-log plot."""
    log_y = np.interp(np.log(x), np.log(xs), np.log(ys))
    return np.exp(log_y)

