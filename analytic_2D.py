#!/Users/cyrus/miniconda3/bin/python3
# from A novel lattice boltzmann model for the poisson equation
from numpy import linspace,array,exp,sinh,cos,pi,fromfunction,zeros,meshgrid
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def read_xml(   f_name,
                fields):
    retDict = {}
    with open(f_name,'r') as f:
        text = f.read().split("\n")
    for line in text:
        for field in fields:
            if "<{}>".format(field[0]) in line:
                temp = list(filter(None,line.split(" ")))
                retDict[field[0]] = field[1](temp[1])
    return(retDict)

def readPalabos(filename,
                nx,
                ny):
    with open(filename,'r') as f:
        text = array(f.read().split(" "))[0:-1].reshape(nx+1,ny+1)
    return(text)

fields = [
            ('lx',int),
            ('ly',int),
            ('resolution',int),
            ('tau_phi',float),
            ('K_0',float),
        ]
params = read_xml("params.xml",fields)

resolution = params['resolution']
nx = int(params['lx']*resolution)
ny = int(params['ly']*resolution)
K_0 = params['K_0']
tau = params['tau_phi']


def u(x_,y_,mu):
    temp = zeros([len(x_),len(y_)])
    for i,x in enumerate(x_):
        for j,y in enumerate(y_):
            temp[i,j] = (cos(pi*x)*sinh(mu*(1-y))/sinh(mu))
    return(temp)

x = linspace(0,1,nx+1)
y = linspace(0,1,ny+1)

lamb = 2
mu = (lamb**2+pi**2)**(1/2)

soln = u(x,y,mu)

lbm_data=readPalabos("concentration_final.dat",nx,ny)
#x2 = linspace(0,1,len(lbm_data))

plt.figure(0, figsize = (6,6))
ax =plt.gca()
ax.contour(x,y,soln.T, colors = 'r', levels = 30 )
ax.contour(x,y,lbm_data.T, colors = 'b', levels = 30 )
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.axis('equal')

yy,xx = meshgrid(x,y)

fig = plt.figure(1, figsize = (12,5))
ax = fig.add_subplot(1,2,1,projection = '3d')
ax.plot_wireframe(xx,yy,soln)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("u(x,y)")

ax = fig.add_subplot(1,2,2,projection = '3d')
ax.plot_wireframe(xx,yy,lbm_data.astype(float), color = 'r')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("u(x,y)")
plt.savefig("figure3.png")
plt.show()
