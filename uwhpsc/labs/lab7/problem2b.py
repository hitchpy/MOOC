
convert_to_cm = False
# Note that in problem2.py,  this was the string 'False' rather than the 
# boolean constant False.   Any nonempty string evaluates to True if
# used in an if test as is done below.

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
t = float(t)
v = velocity(t)

if convert_to_cm:
    print "After %g seconds, the velocity is %g cm/sec"  % (t,v)
else:
    print "After %g seconds, the velocity is %g m/sec"  % (t,v)

