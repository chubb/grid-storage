############
# Flywheel #
############
def PrintMLR(m, l, r, specific_energy):
    specific_energy_Wh_per_kg = specific_energy / c['J_per_kWh'] * 1000.0
    
    print("Cylindrical Flywheel Specs (rotor only):")
    print("mass            = %5.0f kg" % m)
    print("radius          = %5.2f m" % r)
    print("length          = %5.2f m" % l)
    print("specific energy = %5.1f Wh/kg" % specific_energy_Wh_per_kg)
    print("energy density  = %5.1f Wh/L" % (specific_energy_Wh_per_kg * pm_fw['steel']['density'] / 1000.0))

def PrintFlywheelCosts(cost, total_flywheel_icc, rated_energy):
    rated_energy_kWh = rated_energy/c['J_per_kWh']
    
    print("Component          Cost [$]   Cost[$/kWh]")
    print("-----------------------------------------")
    for k in cost:
        print("%-17s %9.2f  %12.2f" % (k, cost[k], cost[k]/rated_energy_kWh))
    print("-----------------------------------------")
    print("TOTAL             %9.2f  %12.2f" % (total_flywheel_icc, total_flywheel_icc/rated_energy_kWh))


############
# Battery  #
############

def PlotBatteryLife60():
    dod = 0.75
    (num_cycles, capacity_at_eol) = BatteryLife60(dod)
    print("""As shown in the plot, when charge cycles are to a depth of discharge
of %2.0f%%, the battery capacity will be 60%% of its initial capacity after %4.0f cycles.""" % 
    (dod*100.0, num_cycles))
    
    dods = np.array([0.01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0])
    (num_cycles, capacity_at_eol) = BatteryLife60(dods)
    #num_cycles = [BatteryLife60(x) for x in dods]
    pylab.semilogy(dods, num_cycles)
    pylab.title("Number of cycles unil End of Life\nWhere End of Life is defined as 60% of initial capacity")
    pylab.grid("on")
    pylab.xlabel("Depth of Discharge [#]")
    pylab.ylabel("Cycles until EoL [#]");


def PlotUpdateBatteryCapacity():
    n = 13000
    capacities = np.zeros(n)
    capacity_0 = 100.0
    capacities[0] = capacity_0
    for i in range(1, n):
        capacities[i] = UpdateBatteryCapacity(capacity_0, capacities[i-1], capacity_0/2.0)

    figure(2)
    pylab.plot(capacities)
    pylab.xlabel("Cycles")
    

################
# Battery Cost #
################
def PlotLiIonCost():
    print("Selected Points from the trend line")
    print("Year    Price [$/kWh]")
    print("%4.0f    %4.0f" % (2017, LiIonCost(2017)*c['J_per_kWh']))
    print("%4.0f    %4.0f" % (2022, LiIonCost(2022)*c['J_per_kWh']))
    x = np.array(range(2005, 2030))
    pylab.figure(figsize=(6,3))
    pylab.plot(x, LiIonCost(x)*c['J_per_kWh'])
    pylab.ylim(0, 600)
    pylab.xlim(2010, 2030)
    pylab.xticks(range(2010, 2030, 2))
    pylab.grid("on")
    pylab.ylabel("Li Ion Cost [$/kWh]")
    pylab.xlabel("Year");