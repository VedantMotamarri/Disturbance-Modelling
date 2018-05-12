import numpy as np
import constants_1U as constants
import qnv as qnv
import satellite as satellite

'''
this function takes input:
    position of centre of mass of satellite in eci frame
    quarternion to convert vector in body frame to eci frame
    inertia matrix of satellite
             gives output:
    torque due to gravity gradient about centre of mass in body frame
'''


def gg_torque(s):
    q = s.getQ()
    qi = qnv.quatInv(q)
    v_pos_sat_i = s.getPos()
    v_pos_sat_b = qnv.quatRotate(qi, v_pos_sat_i)
    pos_norm = np.linalg.norm(v_pos_sat_b)
    m_inertia = constants.m_INERTIA  # moment of inertia matrix in body frame
    m = constants.M
    g = constants.G
    v_t_gg_b = 3 * m * g * (np.cross(v_pos_sat_b, np.matmul(m_inertia, v_pos_sat_b))) / pos_norm ** 5
    return v_t_gg_b


'''
this function takes input:
    velocity of COM of satellite in eci frame
    quarternion to convert a vector in body frame to eci frame
    vector between COM and geometric centre expressed in body frame
              gives output:
    torque due to air drag about COM in body frame
'''


def aero_torque(s):
    q = s.getQ()
    qi = qnv.quatInv(q)
    v_vel_i = s.getVel()
    v_vel_b = qnv.quatRotate(qi, v_vel_i)
    r_com = constants.r_com  # vector containing coordinates of COM of satellite in body frame
    cd = constants.cd  # aerodynamic drag coefficient of satellite
    rho = constants.rho  # density of atmosphere at low earth orbit
    vel_norm = np.linalg.norm(v_vel_b)
    cosw1 = np.dot(v_vel_b, [1, 0, 0]) / x
    cosw2 = np.dot(v_vel_b, [0, 1, 0]) / x
    cosw3 = np.dot(v_vel_b, [0, 0, 1]) / x
    l = constants.length  # length of cube in body frame
    area = l * l * (abs(cosw1) + abs(cosw2) + abs(cosw3))  # area of satellite perpendicular to velocity
    v_t_ad_b = np.cross(r_com, v_vel_b) * rho * cd * vel_norm * area / 2
    return v_t_ad_b


'''
this function takes input:
    sun vector in eci frame
    quarternion to convert a vector in body frame to eci frame
    vector between COM and geometric centre expressed in body frame
              gives output:
    torque due to solar drag about COM in body frame
'''


def solar_torque(s):
    q = s.getQ()
    qi = qnv.quatInv(q)
    r_com = constants.r_com
    e = constants.e  # coefficient of reflectivity of satellite
    p = constants.p  # solar radiation pressure at low earth orbit
    v_sv_i = s.getSun_i()  # unit sun vector in inertial frame obtained from satellite object
    v_sv_b = qnv.quatRotate(qi, v_sv_i) / np.linalg.norm(v_sv_i)
    cosw1 = np.dot(v_sv_b, [1, 0, 0])
    cosw2 = np.dot(v_sv_b, [0, 1, 0])
    cosw3 = np.dot(v_sv_b, [0, 0, 1])
    l = constants.length
    area = l * l * (abs(cosw1) + abs(cosw2) + abs(cosw3))
    # area of satellite perpendicular to sun vector
    v_t_sd1_b = np.cross(r_com, v_sv_b) * p * (1 - e) * area  # torque due to absorption
    v_t_sd2_b = np.cross(r_com, [abs(cosw1) * cosw1, abs(cosw2) * cosw2, abs(cosw3) * cosw3]) * 2 * e * p * l * l
    # reflection torque
    v_t_sd_b = v_t_sd1_b + v_t_sd2_b
    return v_t_sd_b
