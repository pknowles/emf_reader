import serial
import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np
import re
from scipy import linalg
from scipy import optimize

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.grid()

x=[]
y=[]
z=[]

try:
    with open('samples.npy', 'rb') as f:
        data = np.load(f)
        x = list(data[0])
        y = list(data[1])
        z = list(data[2])
except IOError:
    pass



plot = ax.scatter(x, y, z)
plot_cal = ax.scatter(x, y, z)

i=0
try:
    ser = serial.Serial('COM7',9600)
    ser.close()
    ser.open()
except serial.SerialException:
    print("Failed to connect")
    ser = None

from   scipy.spatial import ConvexHull, convex_hull_plot_2d
import numpy as np
from   numpy.linalg import eig, inv
def ls_ellipsoid(xx,yy,zz):                                  
    #finds best fit ellipsoid. Found at http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
    #least squares fit to a 3D-ellipsoid
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz  = 1
    #
    # Note that sometimes it is expressed as a solution to
    #  Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz  = 1
    # where the last six terms have a factor of 2 in them
    # This is in anticipation of forming a matrix with the polynomial coefficients.
    # Those terms with factors of 2 are all off diagonal elements.  These contribute
    # two terms when multiplied out (symmetric) so would need to be divided by two
    
    # change xx from vector of length N to Nx1 matrix so we can use hstack
    x = xx[:,np.newaxis]
    y = yy[:,np.newaxis]
    z = zz[:,np.newaxis]
    
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz = 1
    J = np.hstack((x*x,y*y,z*z,x*y,x*z,y*z, x, y, z))
    K = np.ones_like(x) #column of ones
    
    #np.hstack performs a loop over all samples and creates
    #a row in J for each x,y,z sample:
    # J[ix,0] = x[ix]*x[ix]
    # J[ix,1] = y[ix]*y[ix]
    # etc.
    
    JT=J.transpose()
    JTJ = np.dot(JT,J)
    InvJTJ=np.linalg.inv(JTJ)
    ABC= np.dot(InvJTJ, np.dot(JT,K))

    # Rearrange, move the 1 to the other side
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz - 1 = 0
    #    or
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz + J = 0
    #  where J = -1
    eansa=np.append(ABC,-1)

    return (eansa)

def polyToParams3D(vec,printMe):                             
    #gets 3D parameters of an ellipsoid. Found at http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
    # convert the polynomial form of the 3D-ellipsoid to parameters
    # center, axes, and transformation matrix
    # vec is the vector whose elements are the polynomial
    # coefficients A..J
    # returns (center, axes, rotation matrix)
    
    #Algebraic form: X.T * Amat * X --> polynomial form
    
    if printMe: print('\npolynomial\n',vec)
    
    Amat=np.array(
    [
    [ vec[0],     vec[3]/2.0, vec[4]/2.0, vec[6]/2.0 ],
    [ vec[3]/2.0, vec[1],     vec[5]/2.0, vec[7]/2.0 ],
    [ vec[4]/2.0, vec[5]/2.0, vec[2],     vec[8]/2.0 ],
    [ vec[6]/2.0, vec[7]/2.0, vec[8]/2.0, vec[9]     ]
    ])
    
    if printMe: print('\nAlgebraic form of polynomial\n',Amat)
    
    #See B.Bartoni, Preprint SMU-HEP-10-14 Multi-dimensional Ellipsoidal Fitting
    # equation 20 for the following method for finding the center
    A3=Amat[0:3,0:3]
    A3inv=inv(A3)
    ofs=vec[6:9]/2.0
    center=-np.dot(A3inv,ofs)
    if printMe: print('\nCenter at:',center)
    
    # Center the ellipsoid at the origin
    Tofs=np.eye(4)
    Tofs[3,0:3]=center
    R = np.dot(Tofs,np.dot(Amat,Tofs.T))
    if printMe: print('\nAlgebraic form translated to center\n',R,'\n')
    
    R3=R[0:3,0:3]
    R3test=R3/R3[0,0]
    # print('normed \n',R3test)
    s1=-R[3, 3]
    R3S=R3/s1
    (el,ec)=eig(R3S)
    
    recip=1.0/np.abs(el)
    axes=np.sqrt(recip)
    if printMe: print('\nAxes are\n',axes  ,'\n')
    
    inve=inv(ec) #inverse is actually the transpose here
    if printMe: print('\nRotation matrix\n',inve)
    return (center,axes,inve)

def ellipsoid_fit(data):
    def unpack(x):
        offset=x[:3]
        offset.shape = (3,)
        rotate=x[3:]
        rotate.shape = (3,3)
        return offset, rotate
    def error(x, args):
        offset, rotate = unpack(x)
        s = 0.0
        for point in args.transpose():
            point = np.dot(rotate, point - offset)
            s += (100.0 - point.dot(point)**0.5)**2
        return s
    #x0 = np.array([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1])
    x0 = np.array(
        [32.04076239, 37.97431851, -65.1762087] +
        [2.31681277, 0.21903865, 0.300846858,
            0.21903865, 2.86506899, 0.15303164,
            0.30084685, 0.15303164, 2.47527054]
        )
    result = optimize.minimize(error, x0=x0, args=data)
    print(result)
    return unpack(result.x)

#https://github.com/nliaudat/magnetometer_calibration/blob/main/calibrate.py
class Magnetometer(object):
    '''
        To obtain Gravitation Field (raw format):
    1) get the Total Field for your location from here:
       http://www.ngdc.noaa.gov/geomag-web (tab Magnetic Field)
       es. Total Field = 47,241.3 nT | my val :47'789.7
    2) Convert this values to Gauss (1nT = 10E-5G)
       es. Total Field = 47,241.3 nT = 0.47241G
    3) Convert Total Field to Raw value Total Field, which is the
       Raw Gravitation Field we are searching for
       Read your magnetometer datasheet and find your gain value,
       Which should be the same of the collected raw points
       es. on HMC5883L, given +_ 1.3 Ga as Sensor Field Range settings
           Gain (LSB/Gauss) = 1090 
           Raw Total Field = Gain * Total Field
           0.47241 * 1090 = ~515  |
           
        -----------------------------------------------
         gain (LSB/Gauss) values for HMC5883L
            0.88 Ga => 1370 
            1.3 Ga => 1090 
            1.9 Ga => 820
            2.5 Ga => 660 
            4.0 Ga => 440
            4.7 Ga => 390 
            5.6 Ga => 330
            8.1 Ga => 230 
        -----------------------------------------------
     references :
        -  https://teslabs.com/articles/magnetometer-calibration/      
        -  https://www.best-microcontroller-projects.com/hmc5883l.html
    '''
    MField = 110
    def __init__(self, F=MField): 
        # initialize values
        self.F   = F
        self.b   = np.zeros([3, 1])
        self.A_1 = np.eye(3)
        
    def run(self, data):
        print("shape of data:",data.shape)
        #print("datatype of data:",data.dtype)
        #print("First 5 rows raw:\n", data[:5])
        
        # ellipsoid fit
        M, n, d = self.__ellipsoid_fit(data)

        # calibration parameters
        M_1 = linalg.inv(M)
        self.b = -np.dot(M_1, n)
        self.A_1 = np.real(self.F / np.sqrt(np.dot(n.T, np.dot(M_1, n)) - d) * linalg.sqrtm(M))

        eansa = ls_ellipsoid(data[0],data[1],data[2]) #get ellipsoid polynomial coefficients
        print("coefficients:"  , eansa)
        center,axes,inve = polyToParams3D(eansa,False)
        L = np.diag([100/axes[0],100/axes[1],100/axes[2]])
        #M=np.dot(inve,np.dot(L,linalg.inv(inve)))
        M=np.dot(L.T,inve.T)

        offset, rotate = ellipsoid_fit(data)
        self.b = offset
        self.A_1 = rotate
        
        #print("M:\n", M, "\nn:\n", n, "\nd:\n", d)        
        #print("M_1:\n",M_1, "\nb:\n", self.b, "\nA_1:\n", self.A_1)

        print("Soft iron transformation matrix:\n",self.A_1)
        print("Hard iron bias:\n", self.b)

        result = []
        for row in data.transpose():
        
            # subtract the hard iron offset
            xm_off = row[0]-self.b[0]
            ym_off  = row[1]-self.b[1]
            zm_off  = row[2]-self.b[2]
            
            #multiply by the inverse soft iron offset
            xm_cal = xm_off *  self.A_1[0,0] + ym_off *  self.A_1[0,1]  + zm_off *  self.A_1[0,2] 
            ym_cal = xm_off *  self.A_1[1,0] + ym_off *  self.A_1[1,1]  + zm_off *  self.A_1[1,2] 
            zm_cal = xm_off *  self.A_1[2,0] + ym_off *  self.A_1[2,1]  + zm_off *  self.A_1[2,2] 

            result = np.append(result, np.array([xm_cal, ym_cal, zm_cal]) )#, axis=0 )
            #result = np.append(result, np.array([xm_off, ym_off, zm_off]) )

        result.shape = (-1, 3)
        result = result.transpose()
        return result.reshape(3, -1)

        print("*************************" )        
        print("code to paste : " )
        print("*************************" )  
        print("float hard_iron_bias_x = ", float(self.b[0]), ";")
        print("float hard_iron_bias_y = " , float(self.b[1]), ";")
        print("float hard_iron_bias_z = " , float(self.b[2]), ";")
        print("\n")
        print("double soft_iron_bias_xx = " , float(self.A_1[0,0]), ";")
        print("double soft_iron_bias_xy = " , float(self.A_1[1,0]), ";")
        print("double soft_iron_bias_xz = " , float(self.A_1[2,0]), ";")
        print("\n")
        print("double soft_iron_bias_yx = " , float(self.A_1[0,1]), ";")
        print("double soft_iron_bias_yy = " , float(self.A_1[1,1]), ";")
        print("double soft_iron_bias_yz = " , float(self.A_1[2,1]), ";")
        print("\n")
        print("double soft_iron_bias_zx = " , float(self.A_1[0,2]), ";")
        print("double soft_iron_bias_zy = " , float(self.A_1[1,2]), ";")
        print("double soft_iron_bias_zz = " , float(self.A_1[2,2]), ";")
        print("\n")

    def __ellipsoid_fit(self, s):
        ''' Estimate ellipsoid parameters from a set of points.
            Parameters
            ----------
            s : array_like
              The samples (M,N) where M=3 (x,y,z) and N=number of samples.
            Returns
            -------
            M, n, d : array_like, array_like, float
              The ellipsoid parameters M, n, d.
            References
            ----------
            .. [1] Qingde Li; Griffiths, J.G., "Least squares ellipsoid specific
               fitting," in Geometric Modeling and Processing, 2004.
               Proceedings, vol., no., pp.335-340, 2004
        '''

        # D (samples)
        D = np.array([s[0]**2., s[1]**2., s[2]**2.,
                      2.*s[1]*s[2], 2.*s[0]*s[2], 2.*s[0]*s[1],
                      2.*s[0], 2.*s[1], 2.*s[2], np.ones_like(s[0])])

        # S, S_11, S_12, S_21, S_22 (eq. 11)
        S = np.dot(D, D.T)
        S_11 = S[:6,:6]
        S_12 = S[:6,6:]
        S_21 = S[6:,:6]
        S_22 = S[6:,6:]

        # C (Eq. 8, k=4)
        C = np.array([[-1,  1,  1,  0,  0,  0],
                      [ 1, -1,  1,  0,  0,  0],
                      [ 1,  1, -1,  0,  0,  0],
                      [ 0,  0,  0, -4,  0,  0],
                      [ 0,  0,  0,  0, -4,  0],
                      [ 0,  0,  0,  0,  0, -4]])

        # v_1 (eq. 15, solution)
        E = np.dot(linalg.inv(C),
                   S_11 - np.dot(S_12, np.dot(linalg.inv(S_22), S_21)))

        E_w, E_v = np.linalg.eig(E)

        v_1 = E_v[:, np.argmax(E_w)]
        if v_1[0] < 0: v_1 = -v_1

        # v_2 (eq. 13, solution)
        v_2 = np.dot(np.dot(-np.linalg.inv(S_22), S_21), v_1)

        # quadric-form parameters
        M = np.array([[v_1[0], v_1[3], v_1[4]],
                      [v_1[3], v_1[1], v_1[5]],
                      [v_1[4], v_1[5], v_1[2]]])
        n = np.array([[v_2[0]],
                      [v_2[1]],
                      [v_2[2]]])
        d = v_2[3]

        return M, n, d

mag = Magnetometer()
min_sphere = None
max_sphere = None
first = True

def update(a):
    global first
    updates = 0
    while ser is not None and ser.in_waiting:
        line = ser.readline()
        try:
            line = line.decode().strip(" \n\r")
            if not line.startswith("Raw:"):
                continue
            val = tuple(map(float, line[4:].split(",")))[6:]
        except ValueError:
            continue
        if len(val) != 3:
            continue
        
        if len(x):
            if abs(x[-1] - val[0]) < 2: continue
            if abs(y[-1] - val[1]) < 2: continue
            if abs(z[-1] - val[2]) < 2: continue

        updates += 1

        x.append(val[0])
        y.append(val[1])
        z.append(val[2])

    if not first and updates == 0:
        return ()
    first = False

    limx = (min(x), max(x))
    limy = (min(y), max(y))
    limz = (min(z), max(z))
    #limits = min(min(x), min(y), min(z)), max(max(x), max(y), max(z))
    plot._offsets3d = (x, y, z)

    if len(x) >= 10:
        data = np.array((x, y, z))
        calibrated = mag.run(data)
        #limits = min(limits[0], np.amin(calibrated)), max(limits[0], np.max(calibrated))
        cmin = np.amin(calibrated, axis=1)
        cmax = np.amax(calibrated, axis=1)
        limx = (min(limx[0], cmin[0]), max(limx[1], cmax[0]))
        limy = (min(limy[0], cmin[1]), max(limy[1], cmax[1]))
        limz = (min(limz[0], cmin[2]), max(limz[1], cmax[2]))
        plot_cal._offsets3d = calibrated
    
        lengths = []
        for point in calibrated.transpose():
            lengths += [(point[0]**2 + point[1]**2 + point[2]**2) ** 0.5]
        lengths = sorted(lengths)
        min_r = lengths[len(lengths)//20]
        max_r = lengths[-1-len(lengths)//20]
        print("Fit", min_r, max_r)

        # draw sphere
        def makeSphere(r):
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x = np.cos(u)*np.sin(v)*r
            y = np.sin(u)*np.sin(v)*r
            z = np.cos(v)*r
            return ax.plot_wireframe(x, y, z, color="r")

        global min_sphere, max_sphere
        if min_sphere: min_sphere.remove()
        min_sphere = makeSphere(min_r)
        if max_sphere: max_sphere.remove()
        max_sphere = makeSphere(max_r)


    # Lol. What a fkn mess
    # https://github.com/matplotlib/matplotlib/issues/1077/
    aspect = ax.get_box_aspect()
    aspect = aspect[2] / aspect[0]
    #ax.set_xlim(*limits)
    #ax.set_ylim(*limits)
    #ax.set_zlim(limits[0] * aspect, limits[1] * aspect)
    range_max = max(limx[1] - limx[0], limy[1] - limy[0], limz[1] - limz[0])
    limx = (limx[0], limx[0] + range_max)
    limy = (limy[0], limy[0] + range_max)
    limz = (limz[0], limz[0] + range_max)
    ax.set_xlim(*limx)
    ax.set_ylim(*limy)
    ax.set_zlim(limz[0] * aspect, limz[1] * aspect)

    return (plot,)

ani = matplotlib.animation.FuncAnimation(fig, update, interval=100) 
plt.show()

#with open('samples.npy', 'wb') as f:
#    np.save(f, np.array((x, y, z)))
