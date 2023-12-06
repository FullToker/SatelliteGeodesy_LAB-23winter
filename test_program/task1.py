import numpy as np
import math
from  numpy import  linalg 
import matplotlib.pyplot as plt

# get all the observations(MJD GPSweek and XYZ)
def read_xyz_file(file_path):
    mjd = []
    coors= []

    try:
        with open(file_path, 'r') as file:
            for line in file:
                data = line.strip().split()
                mjd.append((float(data[3]), int(data[4])))
                coors.append((float(data[6]), float(data[7]), float(data[8])))

    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

    return mjd, coors


KIRU_path = 'data/KIRU_ig1.xyz'  
MORP_path = 'data/MORP_ig1.xyz' 
REYK_path = 'data/REYK_ig1.xyz' 
kiru_t, kiru_coors = read_xyz_file(KIRU_path)
morp_t, morp_coors = read_xyz_file(MORP_path)
reyk_t, reyk_coors = read_xyz_file(REYK_path)
# print(kiru_coors[0][0])

# get the position of
def read_snx_file(file_path,sta_name):
    # (Longitude_ _Latitude__ _Height)
    position = []
    long= []
    lat= []

    try:
        with open(file_path, 'r') as file:
            for line in file:
                data = line.strip().split()
                if data and data[0] in sta_name:
                    long=[(float(data[-7]), float(data[-6]), float(data[-5]))]
                    lat=[(float(data[-4]), float(data[-3]), float(data[-2]))]
                    position.append([long, lat, float(data[-1])])

    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

    return position

path1="data/Discontinuities_CONFIRMED.snx"
sta_name=['KIRU', 'MORP', 'REYK']
sta_pos=read_snx_file(path1,sta_name)
print("station LLH:")
print(sta_pos)

# transfer from dms to radian 
def dms_to_ra(degrees, minutes, seconds):
    total_degrees = degrees + minutes / 60 + seconds / 3600
    radians = math.radians(total_degrees)
    
    return radians

# calculate the rotation matrix
def cal_Rota(angle:float,flag:int):
    ro=np.eye(3);
    ra=angle
    if flag==2:
        ro[0][0]=math.cos(ra)
        ro[2][2]=math.cos(ra)
        ro[2][0]=math.sin(ra)
        ro[0][2]=-math.sin(ra)
    elif flag==3:
        ro[0][0]=math.cos(ra)
        ro[1][1]=math.cos(ra)
        ro[1][0]=-math.sin(ra)
        ro[0][1]=math.sin(ra)

    return ro

# transfer time from MJD to year
def mjd2year(_t:list):
    num=len(_t)
    new_t=[]
    for i in range(num):
        jd=_t[i][0]+2400000.5
        T=(jd - 2451545.0) / 36525
        year = 2000 + T * 100
        new_t.append(year)
    return new_t


# read_ITRF2008
def read_ITRF(file_path,sta_name):
    # (Longitude_ _Latitude__ _Height)
    itrf={}

    try:
        with open(file_path, 'r') as file:
            lines=file.readlines()
            i=0
            lines_data=[[],[],[]]
            while(i<len(lines)):
                data = lines[i].strip().split()
                if len(data)==10:
                    if data[-7] == sta_name[0]:
                        lines_data[0].append(lines[i])
                    elif data[-7] == sta_name[1]:
                        lines_data[1].append(lines[i])
                    elif data[-7] == sta_name[2]:
                        lines_data[2].append(lines[i])
                if len(data)>11:
                    if data[-10] == sta_name[0]:
                        lines_data[0].append(lines[i])
                    elif data[-10] == sta_name[1]:
                        lines_data[1].append(lines[i])
                    elif data[-10] == sta_name[2]:
                        lines_data[2].append(lines[i])
                i+=1
        
        i=0
        while(i<3):
            name=sta_name[i]
            lines_ls=lines_data[i]
            con=[]
            for lines in lines_ls:
                data = lines.strip().split()
                if len(data)<11:
                    x=float(data[-6])
                    y=float(data[-5])
                    z=float(data[-4])
                    con.append([{'x':x,'y':y,'z':z}])
                else:
                    x=float(data[-9])
                    y=float(data[-8])
                    z=float(data[-7])
                    start_t=data[-2]
                    end_t=data[-1]
                    con.append([{'x':x,'y':y,'z':z},{'start':start_t,'end':end_t}])
            itrf[name]=con 
            i+=1

    except FileNotFoundError:
        print(f"File not found: {file_path}")
    #except Exception as e:
    #   print(f"An error occurred: {e}")

    return itrf

path1="data/ITRF2008_GNSS.SSC.txt"
# sta_name=['KIRU', 'MORP', 'REYK']
sta_itrf=read_ITRF(path1,sta_name)
print("station ITRF Postion:")
print(sta_itrf)

print(sta_pos[0][0])

# Transform to radian of llh;
sta_llh=[]
for llh in sta_pos:
    mid=[]
    long=dms_to_ra(llh[0][0][0], llh[0][0][1], llh[0][0][2])
    lat=dms_to_ra(llh[1][0][0], llh[1][0][1], llh[1][0][2])
    height=llh[2]
    mid=[long,lat,height]
    sta_llh.append(mid)
print(sta_llh)

print(sta_itrf['KIRU'][0])

# calculate each coors on the local horizontal system
def trans2LHS(sta_itrf:list, sta_llh:list, sta_name:list, flag:int, _t:list,_coors:list):
    name=sta_name[flag]
    # need modify
    mid=sta_itrf[name]
    r2=cal_Rota(-sta_llh[flag][1],2)
    r3=cal_Rota(sta_llh[flag][0],3)
    R=np.dot(r2,r3)
    mid=mid[0]
    x0=mid[0]['x']
    y0=mid[0]['y']
    z0=mid[0]['z']

    pos_LHS=[]
    for coor in _coors:
        delta_m=np.array([[coor[0]-x0], [coor[1]-y0], [coor[2]-z0]])
        new_coor=np.dot(R,delta_m)
        pos_LHS.append(new_coor)
    
    return pos_LHS

coors_lhs_kiru=trans2LHS(sta_itrf,sta_llh,sta_name,0,kiru_t,kiru_coors)
coors_lhs_morp=trans2LHS(sta_itrf,sta_llh,sta_name,1,morp_t,morp_coors)
coors_lhs_reyk=trans2LHS(sta_itrf,sta_llh,sta_name,2,reyk_t,reyk_coors)

print(float(coors_lhs_kiru[0][1]))


#间接平差
def Adjust_example(n,t,A,L,P):
     """
    n 观测值个数；
    t 参数个数；
    A 误差方程系数矩阵数组长度为 n*t，已知；
    L 观测值向量，数组长度为n，已知；
    P 观测值权数组，只有权矩阵的对角线元素数组长度为n，已知
    ZX 参数平差值向量，数组长度为t，待计算
    N 参数平差值的权逆阵Qx; 
    V 观测值残差向量数组长度为n，待计算；
    函数返回值—若计算成功，返回单位权中误差μ。
    """
     N = np.dot(np.dot(A.T,P),A)
     U = np.dot(np.dot(A.T,P),L)
     # print(N)
     # print(U)
     ZX = np.dot(linalg.inv(N),U)
     V = np.dot(A,ZX) - L
     μ = np.sqrt(np.dot(np.dot(V.T,P),V)/(n-t))
     Qx = linalg.inv(N)
     
     return (ZX,V,μ,Qx)


# calculate the parameters of the time series
# flag: 0:north 1:east 2:height
def cal_para(lhs:list,_t:list,flag:int):
    num=len(lhs)
    if len(_t) != num:
        return False
    y=np.empty((num,1))
    A=np.empty((num,4))
    for i in range(num):
        y[i][0]=float(lhs[i][flag])
        A[i][0]=1
        A[i][1]=_t[i]
        A[i][2]=math.cos(2*math.pi*_t[i])
        A[i][3]=math.sin(2*math.pi*_t[i])
    P=np.eye(num)
    print(P[3][3])
    iter=0
    epi=0.00001
    while 1:
        zx, v, _, _=Adjust_example(num,4,A,y,P)
        iter+=1
        for i in range(len(y)):
            y[i]+=v[i]
        if abs(np.max(v))<epi and abs(np.min(v))<epi:
            break
    print(iter)
    return y,zx


new_kiru_t=mjd2year(kiru_t)
new_lhs_kiru_north, kiru_Nbeta=cal_para(coors_lhs_kiru,new_kiru_t,0)
new_lhs_kiru_east, kiru_Ebeta=cal_para(coors_lhs_kiru,new_kiru_t,1)
new_lhs_kiru_h, kiru_Hbeta=cal_para(coors_lhs_kiru,new_kiru_t,2)

print("parameter is :")
print(kiru_Nbeta)


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(8, 10))
ax1.plot(new_kiru_t, new_lhs_kiru_north)
ax1.set_ylabel('north')
ax1.legend()

ax2.plot(new_kiru_t, new_lhs_kiru_east)
ax2.set_ylabel('east')
ax2.legend()

ax3.plot(new_kiru_t, new_lhs_kiru_h)
ax3.set_ylabel('height')
ax3.legend()

fig.suptitle('Three Subplots with Shared x-axis')
plt.tight_layout()
plt.show()
'''

delta_kiru_E=[]
for i in range(len(new_lhs_kiru_north)):
    #mid=new_lhs_kiru_east[i]-kiru_Nbeta[1]*new_kiru_t[i]
    #mid=new_lhs_kiru_north[i]-kiru_Nbeta[0]
    mid=kiru_Nbeta[1]*new_kiru_t[i]
    delta_kiru_E.append(mid)

plt.plot(new_kiru_t, delta_kiru_E)
plt.xlabel('t')
plt.ylabel('north')
plt.title('north')
plt.legend()
plt.show()
'''

    




