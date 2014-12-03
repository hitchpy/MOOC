
convert_to_cm = 'False'

g = 9.81  # gravitational constant, meters/sec**2

if convert_to_cm:
    g = 100.*g  # convert to cm/sec**2

def velocity(t):
    """ 
    Returns velocity of a falling object after t seconds if starting from rest
    """
    v = g*t
    return v

t = raw_input("Input time t = ")
v = velocity(t)

if convert_to_cm:
    print "After %g seconds, the velocity is %g cm/sec"  % v
else:
    print "After %g seconds, the velocity is %g m/sec"  % v

