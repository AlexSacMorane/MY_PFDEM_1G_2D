#------------------------------------------------------------------------------------------------------------------------------------------ #
# Librairies
#------------------------------------------------------------------------------------------------------------------------------------------#

import numpy as np

#------------------------------------------------------------------------------------------------------------------------------------------ #
# Parameters
#------------------------------------------------------------------------------------------------------------------------------------------#

def get_parameters():
    '''
    Define the parameters used in the simulation.
    '''
    #---------------------------------------------------------------------#
    # Norrmalization
    n_dist = 100*1e-6 # m
    n_time = 24*60*60 # s
    n_mol = 0.73*1e3 * n_dist**3 # mol

    #---------------------------------------------------------------------#
    # PFDEM

    n_DEMPF_ite = 150 # number of PFDEM iterations
    n_proc = 4 # number of processors used
    j_total = 0 # index global of results
    save_simulation = False # indicate if the simulation is saved
    n_max_vtk_files = 10 # maximum number of vtk files (can be None to save all files)

    # Select Figures to plot
    # Available:
    # contact_pressure, contact_distrib_m_ed, contact_point_ed, contact_volume, contact_nb_node
    # contact_detection, contact_h_s_v, contact_dem, saturation, m_c_pore
    # m_ed, dt_PF, IC, processor, sphericities, force_applied, maps, n_grain_kc_map, dim_dom
    # shape_evolution, n_vertices, sum_etai_c, mean_etai_c, mass_loss, performances
    # disp_strain, disp_strain_andrade, sample_height, y_contactPoint
    L_figures = ['maps', 'shape_evolution', 'contact_pressure', 'saturation', 'm_c_pore', 'disp_strain']

    # Figure (plot all or current)
    # The maps configuration
    print_all_map_config = False # else only the current one is printed

    #---------------------------------------------------------------------#
    # DEM (Yade)

    # steady state detection
    n_ite_max = 5000 # maximum number of iteration during a DEM step
    n_steady_state_detection = 100 # number of iterations considered in the window
    # the difference between max and min < tolerance * force_applied
    steady_state_detection = 0.02
    # + the force applied must be contained in this window

    # sollicitation
    force_applied = 100e6  # (kg m s-2)
    control_force = True # Boolean to determine if the force is controled with the contact volume

    # DEM material parameters
    # Young modulus
    E = 2e9 # (kg m-1 s-2)

    # Poisson ratio
    Poisson = 0.3

    # Figure (plot all or current)
    # The evolution of the overlap durig DEM steps
    print_all_contact_dem = False # else only the current one is printed
    # The evolution of shape (compared to the initial state)
    print_all_shape_evolution = False # else only the current one is printed

    #---------------------------------------------------------------------#
    # Grain description

    # shape of the grain
    # Sphere, Hex or Proxy_Hex
    Shape = 'Sphere'

    # the radius of grains
    radius = 100*1e-6/n_dist # m/m
    # discretization of the grain
    n_phi = 30

    # size of the indenter
    plate = 0.1*radius # m/m

    #---------------------------------------------------------------------#
    # Phase-Field (Moose)

    # determine if remeshing is available
    remesh = True
    # mesh
    check_database = True
    # remesh
    if remesh:
        size_x_mesh = radius/40
        size_y_mesh = size_x_mesh
        m_size_mesh = (size_x_mesh+size_y_mesh)/2
        margin_mesh_domain = 15*size_x_mesh
    # constant mesh
    else :
        x_min = -1.3*radius
        x_max =  1.3*radius
        y_min = -2.3*radius
        y_max =  2.3*radius
        n_mesh_x = 100
        n_mesh_y = 200
        m_size_mesh = ((x_max-x_min)/(n_mesh_x-1)+(y_max-y_min)/(n_mesh_y-1))/2

    # PF material parameters
    # the energy barrier
    Energy_barrier = 1
    # number of mesh in the interface
    n_int = 6
    # the interface thickness
    w = m_size_mesh*n_int
    # the gradient coefficient
    kappa_eta = Energy_barrier*w*w/9.86
    # the mobility
    Mobility_eff = 3*(100*1e-6/(24*60*60))/(n_dist/n_time) # m.s-1/(m.s-1)

    # temperature
    temperature = 623 # K 
    # molar volume
    V_m = (2.2*1e-5)/(n_dist**3/n_mol) # (m3 mol-1)/(m3 mol-1)
    # constant
    R_cst = (8.32)/(n_dist**2/(n_time**2*n_mol)) # (kg m2 s-2 mol-1 K-1)/(m2 s-2 mol-1)

    # kinetics of dissolution and precipitation
    # it affects the tilting coefficient in Ed
    k_diss = 1*(0.005)/(m_size_mesh) # ed_j = ed_i*m_i/m_j
    k_prec = k_diss*16 # -

    # molar concentration at the equilibrium
    C_eq = (0.73*1e3)/(n_mol/n_dist**3) # (mol m-3)/(mol m-3)

    # diffusion of the solute
    size_film = m_size_mesh*5
    D_solute = (4e-14/2/size_film)/(n_dist*n_dist/n_time) # (m2 s-1)/(m2 s-1)
    n_struct_element = int(round(size_film/m_size_mesh,0))
    struct_element = np.array(np.ones((n_struct_element,n_struct_element)), dtype=bool) # for dilation

    # Aitken method
    # the time stepping and duration of one PF simualtion
    # level 0
    dt_PF_0 = (0.01*24*60*60)/n_time # time step
    # level 1
    dt_PF_1 = dt_PF_0
    m_ed_contact_1 = 0
    # level 2
    dt_PF_2 = dt_PF_1
    m_ed_contact_2 = 0
    # n_t_PF*dt_PF gives the total time duration
    n_t_PF = 200 # number of iterations

    # the criteria on residual
    crit_res = 1e-3
    
    # Contact box detection
    eta_contact_box_detection = 0.1 # value of the phase field searched to determine the contact box

    # Figure (plot all or current)
    # The detection of the contact by a box
    print_all_contact_detection = False # else only the current one is printed

    #---------------------------------------------------------------------#
    # trackers

    L_displacement = []
    L_sum_eta_1 = []
    L_sum_eta_2 = []
    L_sum_c = []
    L_sum_mass = []
    L_m_eta_1 = []
    L_m_eta_2 = []
    L_m_c = []
    L_m_mass = []
    L_distance_extrema = []
    L_equivalent_area = []
    L_contact_overlap = []
    L_contact_area = []
    L_contact_volume_yade = []
    L_contact_volume_moose = []
    L_contact_volume_box = []
    L_t_pf_to_dem_1 = []
    L_t_pf_to_dem_2 = []
    L_t_dem = []
    L_t_dem_to_pf = []
    L_t_pf = []
    L_P_applied = []
    L_n_v_1 = []
    L_n_v_2 = []
    L_n_v_1_target = []
    L_n_v_2_target = []
    L_m_ed = []
    L_m_ed_contact = []
    L_m_ed_large_contact = []
    L_m_ed_plus_contact = []
    L_m_ed_minus_contact = []
    L_m_ed_plus_large_contact = []
    L_m_ed_minus_large_contact = []
    L_ed_contact_point = []
    L_vertices_1_init = None
    L_vertices_1_006 = None
    L_force_applied = []
    L_dt_PF = []
    L_AreaSphericity = []
    L_DiameterSphericity = []
    L_CircleRatioSphericity = []
    L_PerimeterSphericity = []
    L_WidthToLengthRatioSpericity = []
    L_grain_kc_map = []
    L_sample_height = []
    L_y_contactPoint = []
    L_loss_move_pf_eta1 = []
    L_loss_move_pf_c = []
    L_loss_move_pf_m = []
    L_loss_kc_eta1 = []
    L_loss_kc_c = []
    L_loss_kc_m = []
    L_loss_pf_eta1 = []
    L_loss_pf_c = []
    L_loss_pf_m = []
    L_L_profile_sat = []
    L_L_x = []
    L_m_c_pore = []
    L_m_sat_contact = []
    if remesh:
        L_x_min_dom = []
        L_x_max_dom = []
        L_y_min_dom = []
        L_y_max_dom = []
        L_delta_x_max = []
        L_delta_y_max = []

    #---------------------------------------------------------------------#
    # dictionnary

    dict_user = {
    'n_dist': n_dist,
    'n_time': n_time,
    'n_mol': n_mol,
    'n_DEMPF_ite': n_DEMPF_ite,
    'n_proc': n_proc,
    'j_total': j_total,
    'save_simulation': save_simulation,
    'n_max_vtk_files': n_max_vtk_files,
    'L_figures': L_figures,
    'print_all_map_config': print_all_map_config,
    'n_ite_max': n_ite_max,
    'n_steady_state_detection': n_steady_state_detection,
    'steady_state_detection': steady_state_detection,
    'force_applied': force_applied,
    'force_applied_target': force_applied,
    'control_force': control_force,
    'E': E,
    'Poisson': Poisson,
    'print_all_contact_dem': print_all_contact_dem,
    'print_all_shape_evolution': print_all_shape_evolution,
    'Shape': Shape,
    'radius': radius,
    'n_phi': n_phi,
    'plate': plate,
    'remesh': remesh,
    'check_database': check_database,
    'Energy_barrier': Energy_barrier,
    'n_int': n_int,
    'w_int': w,
    'kappa_eta': kappa_eta,
    'Mobility_eff': Mobility_eff,
    'temperature': temperature,
    'V_m': V_m,
    'R_cst': R_cst,
    'k_diss': k_diss,
    'k_prec': k_prec,
    'C_eq': C_eq,
    'size_film': size_film,
    'D_solute': D_solute,
    'struct_element': struct_element,
    'dt_PF_0': dt_PF_0,
    'dt_PF_1': dt_PF_1,
    'Aitken_1': m_ed_contact_1,
    'dt_PF_2': dt_PF_2,
    'Aitken_2': m_ed_contact_2,
    'n_t_PF': n_t_PF,
    'crit_res': crit_res,
    'eta_contact_box_detection': eta_contact_box_detection,
    'print_all_contact_detection': print_all_contact_detection,
    'L_displacement': L_displacement,
    'L_sum_eta_1': L_sum_eta_1,
    'L_sum_eta_2': L_sum_eta_2,
    'L_sum_c': L_sum_c,
    'L_sum_mass': L_sum_mass,
    'L_m_eta_1': L_m_eta_1,
    'L_m_eta_2': L_m_eta_2,
    'L_m_c': L_m_c,
    'L_m_mass': L_m_mass,
    'L_distance_extrema': L_distance_extrema,
    'L_equivalent_area': L_equivalent_area,
    'L_contact_overlap': L_contact_overlap,
    'L_contact_area': L_contact_area,
    'L_contact_volume_yade': L_contact_volume_yade,
    'L_contact_volume_moose': L_contact_volume_moose,
    'L_contact_volume_box': L_contact_volume_box,
    'L_t_pf_to_dem_1': L_t_pf_to_dem_1,
    'L_t_pf_to_dem_2': L_t_pf_to_dem_2,
    'L_t_dem': L_t_dem,
    'L_t_dem_to_pf': L_t_dem_to_pf,
    'L_t_pf': L_t_pf,
    'L_P_applied': L_P_applied,
    'L_n_v_1': L_n_v_1,
    'L_n_v_2': L_n_v_2,
    'L_n_v_1_target': L_n_v_1_target,
    'L_n_v_2_target': L_n_v_2_target,
    'L_m_ed': L_m_ed,
    'L_m_ed_contact': L_m_ed_contact,
    'L_m_ed_large_contact': L_m_ed_large_contact,
    'L_m_ed_plus_contact': L_m_ed_plus_contact,
    'L_m_ed_minus_contact': L_m_ed_minus_contact,
    'L_m_ed_plus_large_contact': L_m_ed_plus_large_contact,
    'L_m_ed_minus_large_contact': L_m_ed_minus_large_contact,
    'L_ed_contact_point': L_ed_contact_point,
    'L_vertices_1_init': L_vertices_1_init,
    'L_vertices_1_006': L_vertices_1_006,
    'L_force_applied': L_force_applied,
    'L_dt_PF': L_dt_PF,
    'L_AreaSphericity': L_AreaSphericity,
    'L_DiameterSphericity': L_DiameterSphericity,
    'L_CircleRatioSphericity': L_CircleRatioSphericity,
    'L_PerimeterSphericity': L_PerimeterSphericity,
    'L_WidthToLengthRatioSpericity': L_WidthToLengthRatioSpericity,
    'L_grain_kc_map': L_grain_kc_map,
    'L_sample_height': L_sample_height,
    'L_y_contactPoint': L_y_contactPoint,
    'L_loss_move_pf_eta1': L_loss_move_pf_eta1,
    'L_loss_move_pf_c': L_loss_move_pf_c,
    'L_loss_move_pf_m': L_loss_move_pf_m,
    'L_loss_kc_eta1': L_loss_kc_eta1,
    'L_loss_kc_c': L_loss_kc_c,
    'L_loss_kc_m': L_loss_kc_m,
    'L_loss_pf_eta1': L_loss_pf_eta1,
    'L_loss_pf_c': L_loss_pf_c,
    'L_loss_pf_m': L_loss_pf_m,
    'L_L_profile_sat': L_L_profile_sat,
    'L_L_x': L_L_x,
    'L_m_c_pore': L_m_c_pore,
    'L_m_sat_contact': L_m_sat_contact
    }

    # specific inputs
    if remesh:
        dict_user['size_x_mesh'] = size_x_mesh
        dict_user['size_y_mesh'] = size_y_mesh
        dict_user['margin_mesh_domain'] = margin_mesh_domain
        dict_user['L_x_min_dom'] = L_x_min_dom
        dict_user['L_x_max_dom'] = L_x_max_dom
        dict_user['L_y_min_dom'] = L_y_min_dom
        dict_user['L_y_max_dom'] = L_y_max_dom
        dict_user['L_delta_x_max'] = L_delta_x_max
        dict_user['L_delta_y_max'] = L_delta_y_max
    else :
        dict_user['x_min'] = x_min
        dict_user['x_max'] = x_max
        dict_user['y_min'] = y_min
        dict_user['y_max'] = y_max
        dict_user['n_mesh_x'] = n_mesh_x
        dict_user['n_mesh_y'] = n_mesh_y

    return dict_user
