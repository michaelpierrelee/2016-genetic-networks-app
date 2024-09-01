# coding: utf8

from . import *
    
def progress(i, n, m, e):
    """
    display a progress bar in the I/O interfacce
    i = degree of progress
    n = maximal value of i
    e = % of progress
    """
    d = ""
    x = 100. / (n + 1)
    
    m += x
    if m >= 1:
        for j in range(int(m)): e += 1
        m -= int(m)
        
    if i == n and round(m,0) >=1 and e < 100: e += 1

    now = datetime.datetime.now()
    d = "\r>Progress: "+ str(e) + "% - " + now.strftime('%H:%M:%S')
    if i == n: d += " ---End---\n"
    try: sys.stdout.write(d)
    except:
        time.sleep(0.001)
        sys.stdout.write(d)
    sys.stdout.flush()
    return m, e
    
