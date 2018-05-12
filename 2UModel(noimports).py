def gg_t(s):


q=s.getQ()
qi=qnv.quatInv(q)
v_pos_sat_i=s.getPos()
v_pos_sat_b=qnv.quatRotate(qi,v_pos_sat_i)
x=np.linalg.norm(v_pos_sat_b)
y=s.m_INERTIA     #moment of inertia matrix in body frame
v_gg_t_b=3*(constants.M)*(constants.G)*(np.cross(v_pos_sat_b,np.matmul(y,v_pos_sat_b)))/(x**5)
return v_gg_t_b


def ad_t(s):


q=s.getQ()
qi=qnv.quatInv(q)
v_vel_sat_i=s.getVel()
v_vel_sat_b=qnv.quatRotate(qi,v_vel_sat_i)
r_com=s.r_com      #vector containing coordinates of COM of satellite in body frame
cd=s.cd            #aerodynamic drag coefficient of satellite
rho=constants.rho  #density of atmosphere at low earth orbit
x=np.linalg.norm(v_vel_sat_b)
cosw1=np.dot(v_vel_sat_b,[1,0,0])/x
cosw2=np.dot(v_vel_sat_b,[0,1,0])/x
cosw3=np.dot(v_vel_sat_b,[0,0,1])/x
l1=s.l1            #length of cuboidal satellite in x direction in body frame
l2=s.l2            #length of cuboidal satellite in y direction in body frame
l3=s.l3            #length of cuboidal satellite in z direction in body frame
area=l2*l3*abs(cosw1)+l1*l3*abs(cosw2)+l1*l2*abs(cosw3) #area of satellite perpendicular to velocity
v_ad_t_b=np.cross(r_com,v_vel_sat_b)*rho*cd*x*area/2
return v_ad_t_b


def sd_t(s):


q=s.getQ()
qi=qnv.quatInv(q)
r_com=s.r_com
e=s.e              #coefficient of reflectivity of satellite
p=constants.p      #solar radiation pressure at low earth orbit
v_sun_i=s.getSun() #unit sun vector in inertial frame obtained from satellite object
v_sun_b=qnv.quatRotate(qi,v_sun_i)/np.linalg.norm(v_sun_i)
cosw1=np.dot(v_sun_b,[1,0,0])
cosw2=np.dot(v_sun_b,[0,1,0])
cosw3=np.dot(v_sun_b,[0,0,1])
l1=s.l1            #length of cuboidal satellite in x direction in body frame
l2=s.l2            #length of cuboidal satellite in y direction in body frame
l3=s.l3            #length of cuboidal satellite in z direction in body frame
area=l2*l3*abs(cosw1)+l1*l3*abs(cosw2)+l1*l2*abs(cosw3) #area of satellite perpendicular to sun vector
v_sd_t_b_abs = np.cross(r_com,v_sun_b)*p*(1-e)*area     #torque due to absorption
v_sd_t_b_ref = np.cross(r_com,[l2*l3*abs(cosw1)*cosw1,l1*l3*abs(cosw2)*cosw2,l1*l2*abs(cosw3)*cosw3])*2*e*p # reflection torque
v_sd_t_b=v_sd_t_b_abs + v_sd_t_b_ref
return v_sd_t 

