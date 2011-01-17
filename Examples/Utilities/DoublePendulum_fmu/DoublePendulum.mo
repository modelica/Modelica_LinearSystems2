within Modelica_LinearSystems2.Examples.Utilities.DoublePendulum_fmu;
model DoublePendulum

  //	Interfaces
type Modelica_Blocks_Interfaces_RealOutput = Real;
type Modelica_Blocks_Interfaces_RealInput = Real;

  //	arbitrary parameters
type Modelica_Mechanics_MultiBody_Types_RotationSequence = Integer(max=3, min=1);
type Modelica_Mechanics_MultiBody_Types_Color = Integer(max=255, min=0);
type Modelica_Mechanics_MultiBody_Types_SpecularCoefficient = Real;
type Modelica_Mechanics_MultiBody_Frames_Quaternions_Orientation = Real;
type Modelica_Mechanics_MultiBody_Types_ShapeExtra = Real;
type Modelica_Mechanics_MultiBody_Types_Axis = Real(unit="1");
type Modelica_Mechanics_MultiBody_Types_AxisLabel = String;
type Modelica_Mechanics_MultiBody_Types_ShapeType = String;

  //	Unit types
type Modelica_SIunits_Distance = Real(min=0.0, unit="m", quantity="Length");
type Modelica_SIunits_Torque = Real(unit="N.m", quantity="Torque");
type Modelica_SIunits_Mass = Real(min=0.0, unit="kg", quantity="Mass");
type Modelica_SIunits_AngularAcceleration = Real(unit="rad/s2", quantity="AngularAcceleration");
type Modelica_SIunits_RotationalDampingConstant = Real(unit="N.m.s/rad", quantity="RotationalDampingConstant");
type Modelica_SIunits_Inertia = Real(unit="kg.m2", quantity="MomentOfInertia");
type Modelica_SIunits_Acceleration = Real(unit="m/s2", quantity="Acceleration");
type Modelica_SIunits_Position = Real(unit="m", quantity="Length");
type Modelica_SIunits_AngularVelocity = Real(unit="rad/s", quantity="AngularVelocity");
type Modelica_SIunits_Power = Real(unit="W", quantity="Power");
type Modelica_SIunits_Angle = Real(displayUnit="deg", unit="rad", quantity="Angle");
type Modelica_SIunits_Velocity = Real(unit="m/s", quantity="Velocity");
type Modelica_SIunits_Force = Real(unit="N", quantity="Force");
type Modelica_SIunits_Density = Real(min=0.0, unit="kg/m3", quantity="Density");
type Modelica_SIunits_TranslationalDampingConstant = Real(unit="N.s/m", quantity="TranslationalDampingConstant");
type Modelica_SIunits_Diameter = Real(min=0.0, unit="m", quantity="Length");
type Modelica_SIunits_Length = Real(unit="m", quantity="Length");

  //Enumerations
type Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB = enumeration(
      world "Resolve in world frame",
      frame_a "Resolve in frame_a",
      frame_b "Resolve in frame_b",
      frame_resolve
        "Resolve in frame_resolve (frame_resolve must be connected)");
type Modelica_Mechanics_MultiBody_Types_GravityTypes = enumeration(
      NoGravity "No gravity field",
      UniformGravity "Uniform gravity field",
      PointGravity "Point gravity field");
type StateSelect = enumeration(
      never,
      avoid,
      default,
      prefer,
      always);

    Modelica.Blocks.Interfaces.RealInput u(start=_u_start);

    Modelica.Blocks.Interfaces.RealOutput w1(unit="1/s");
    Modelica.Blocks.Interfaces.RealOutput phi1;
    Modelica.Blocks.Interfaces.RealOutput s;
    Modelica.Blocks.Interfaces.RealOutput v(unit="m/s");
    Modelica.Blocks.Interfaces.RealOutput w(unit="1/s");
    Modelica.Blocks.Interfaces.RealOutput phi;

    parameter Real _u_start = 0.0;
    parameter Real _rev_phi_start = -1.39626340159546
    "Relative rotation angle from frame_a to frame_b";
    parameter Real _prismatic_v_start = 0
    "First derivative of s (relative velocity)";
    parameter Real _prismatic_s_start = 0
    "Relative distance between frame_a and frame_b";

    parameter Modelica_SIunits_Inertia bodyShape_I_31(min=-1E+060) = 0
    " (3,1) element of inertia tensor";
    parameter Modelica_SIunits_Inertia bodyShape_I_33(min=0.0) = 0.001
    " (3,3) element of inertia tensor";
    parameter Real bodyShape_w_0_start_2_(unit="rad/s") = 0
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real bodyCylinder1_z_0_start_1_(unit="rad/s2") = 0
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real add_k1 = 1 "Gain of upper input";
    parameter Modelica_SIunits_Distance bodyCylinder_innerDiameter = 0
    "Inner diameter of cylinder (0 <= innerDiameter <= Diameter)";
    parameter Real body_r_CM_3_(unit="m") = 0
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Real bodyCylinder1_w_0_start_3_(unit="rad/s") = 0
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Modelica_SIunits_Angle rev_fixed_phi0 = 0
    "Fixed offset angle of housing";
    parameter Real bodyShape_w_0_start_3_(unit="rad/s") = 0
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real body_r_CM_2_(unit="m") = 0
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Real bodyCylinder1_angles_start_1_(displayUnit="deg", unit="rad") = 0
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Real body_w_0_start_3_(unit="rad/s") = 0
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Modelica_SIunits_Density bodyCylinder_density = 900
    "Density of cylinder (e.g., steel: 7700 .. 7900, wood : 400 .. 800)";
    parameter Real bodyCylinder1_z_0_start_3_(unit="rad/s2") = 0
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real bodyShape_r_CM_2_(unit="m") = 0
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Real body_w_0_start_2_(unit="rad/s") = 0
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real body_r_CM_1_(unit="m") = 0
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Real bodyCylinder1_z_0_start_2_(unit="rad/s2") = 0
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real bodyCylinder_angles_start_2_(displayUnit="deg", unit="rad") = 0
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Modelica_SIunits_TranslationalDampingConstant damper1_d(min=0.0) = 0
    "Damping constant";
    parameter Real world_defaultWidthFraction = 20
    "Default for shape width as a fraction of shape length (e.g., for Parts.FixedTranslation)";
    parameter Modelica_SIunits_Inertia bodyShape_I_22(min=0.0) = 0.001
    " (2,2) element of inertia tensor";
    parameter Modelica_SIunits_Inertia bodyShape_I_21(min=-1E+060) = 0
    " (2,1) element of inertia tensor";
    parameter Real body_w_0_start_1_(unit="rad/s") = 0
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real world_mue(unit="m3/s2", min=0.0) = 398600000000000.0
    "Gravity field constant (default = field constant of earth)";
    parameter Real body_angles_start_3_(displayUnit="deg", unit="rad") = 0
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Real bodyCylinder1_w_0_start_1_(unit="rad/s") = 0
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real world_defaultSpecularCoefficient(min=0.0) = 0.7
    "Default reflection of ambient light (= 0: light is completely absorbed)";
    parameter Real bodyShape_z_0_start_1_(unit="rad/s2") = 0
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Modelica_Mechanics_MultiBody_Types_ShapeExtra bodyShape_frameTranslation_extra = 0.0
    " Additional parameter depending on shapeType (see docu of Visualizers.Advanced.Shape).";
    parameter Real world_defaultNm_to_m(unit="N.m/m", min=0.0) = 1000
    "Default scaling of torque arrows (length = torque/defaultNm_to_m)";
    parameter Modelica_SIunits_Inertia bodyShape_I_32(min=-1E+060) = 0
    " (3,2) element of inertia tensor";
    parameter Real bodyCylinder_z_0_start_2_(unit="rad/s2") = 0
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Modelica_SIunits_RotationalDampingConstant damper_d(min=0.0) = 0
    "Damping constant";
    parameter Modelica_Mechanics_MultiBody_Types_ShapeExtra bodyShape_extra = 0.0
    " Additional parameter depending on shapeType (see docu of Visualizers.Advanced.Shape).";
    parameter Real body_angles_start_2_(displayUnit="deg", unit="rad") = 0
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Real bodyCylinder1_w_0_start_2_(unit="rad/s") = 0
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real world_defaultN_to_m(unit="N/m", min=0.0) = 1000
    "Default scaling of force arrows (length = force/defaultN_to_m)";
    parameter Real bodyCylinder_z_0_start_3_(unit="rad/s2") = 0
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real body_angles_start_1_(displayUnit="deg", unit="rad") = 0
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Real bodyCylinder_w_0_start_3_(unit="rad/s") = 0
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real bodyShape_z_0_start_3_(unit="rad/s2") = 0
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Boolean world_driveTrainMechanics3D = false
    "= true, if 3-dim. mechanical effects of Parts.Mounting1D/Rotor1D/BevelGear1D shall be taken into account";
    parameter Real bodyShape_angles_start_2_(displayUnit="deg", unit="rad") = 0
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Modelica_SIunits_Distance bodyCylinder1_diameter = 0.05
    "Diameter of cylinder";
    parameter Modelica_SIunits_Inertia bodyShape_I_11(min=0.0) = 0.001
    " (1,1) element of inertia tensor";
    parameter Real bodyShape_z_0_start_2_(unit="rad/s2") = 0
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real add_k2 = 1 "Gain of lower input";
    parameter Real world_gravityArrowTail_1_(unit="m") = 0
    "Position vector from origin of world frame to arrow tail, resolved in world frame";
    parameter Real bodyCylinder_w_0_start_2_(unit="rad/s") = 0
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real bodyShape_angles_start_3_(displayUnit="deg", unit="rad") = 0
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Modelica_SIunits_Mass m_trolley = 5;
    parameter Modelica_SIunits_Mass m_load = 20;
    parameter Real body_z_0_start_3_(unit="rad/s2") = 0
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real bodyShape_r_CM_3_(unit="m") = 0
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Modelica_SIunits_AngularVelocity w1_start = 0.0;
    parameter Real add1_k2 = 1 "Gain of lower input";
    parameter Real add1_k1 = 1 "Gain of upper input";
    parameter Real world_gravityArrowTail_3_(unit="m") = 0
    "Position vector from origin of world frame to arrow tail, resolved in world frame";
    parameter Real world_defaultFrameDiameterFraction = 40
    "Default for arrow diameter of a coordinate system as a fraction of axis length";
    parameter Real body_z_0_start_2_(unit="rad/s2") = 0
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real bodyCylinder_z_0_start_1_(unit="rad/s2") = 0
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real bodyShape_angles_start_1_(displayUnit="deg", unit="rad") = 0
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Modelica_SIunits_Inertia body_I_11(min=0.0) = 0.001
    " (1,1) element of inertia tensor";
    parameter Real bodyCylinder_angles_start_1_(displayUnit="deg", unit="rad") = 0
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Modelica_SIunits_Angle relativeAngles1_guessAngle1 = 0
    "Select angles[1] such that abs(angles[1] - guessAngle1) is a minimum";
    parameter Real bodyCylinder_w_0_start_1_(unit="rad/s") = 0
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Modelica_SIunits_AngularVelocity w2_start = 0.0;
    parameter Modelica_SIunits_Distance bodyCylinder_diameter = 0.05
    "Diameter of cylinder";
    parameter Real world_gravityArrowTail_2_(unit="m") = 0
    "Position vector from origin of world frame to arrow tail, resolved in world frame";
    parameter Real const_k = 1.5707963267949 "Constant output value";
    parameter Modelica_SIunits_Acceleration world_g = 9.81
    "Constant gravity acceleration";
    parameter Real bodyShape_r_CM_1_(unit="m") = 0
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Modelica_SIunits_Angle revolute2_fixed_phi0 = 0
    "Fixed offset angle of housing";
    parameter Modelica_SIunits_Density bodyCylinder1_density = 900
    "Density of cylinder (e.g., steel: 7700 .. 7900, wood : 400 .. 800)";
    parameter Modelica_SIunits_Inertia body_I_33(min=0.0) = 0.001
    " (3,3) element of inertia tensor";
    parameter Modelica_SIunits_Inertia body_I_32(min=-1E+060) = 0
    " (3,2) element of inertia tensor";
    parameter Modelica_SIunits_Inertia body_I_31(min=-1E+060) = 0
    " (3,1) element of inertia tensor";
    parameter Real bodyCylinder_angles_start_3_(displayUnit="deg", unit="rad") = 0
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Real body_z_0_start_1_(unit="rad/s2") = 0
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Modelica_SIunits_Angle relativeAngles_guessAngle1 = 0
    "Select angles[1] such that abs(angles[1] - guessAngle1) is a minimum";
    parameter Modelica_SIunits_Position prismatic_fixed_s0 = 0
    "Fixed offset position of housing";
    parameter Modelica_SIunits_Diameter world_gravitySphereDiameter = 12742000
    "Diameter of sphere representing gravity center (default = mean diameter of earth)";
    parameter Real bodyCylinder1_angles_start_2_(displayUnit="deg", unit="rad") = 0
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Modelica_SIunits_Inertia body_I_21(min=-1E+060) = 0
    " (2,1) element of inertia tensor";
    parameter Modelica_SIunits_Inertia body_I_22(min=0.0) = 0.001
    " (2,2) element of inertia tensor";
    parameter Modelica_SIunits_Distance bodyCylinder1_innerDiameter = 0
    "Inner diameter of cylinder (0 <= innerDiameter <= Diameter)";
    parameter Modelica_SIunits_Torque revolute2_constantTorque_tau_constant = 0
    "Constant torque (if negative, torque is acting as load)";
    parameter Real bodyShape_w_0_start_1_(unit="rad/s") = 0
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real const1_k = 0 "Constant output value";
    parameter Real bodyCylinder1_angles_start_3_(displayUnit="deg", unit="rad") = 0
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Modelica_SIunits_Angle phi2_start = 10;

  // Parameters with default values
    parameter Real relativePosition_relativePosition_frame_b_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder1_body_R_start_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyCylinder_body_cylinderColor_1_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real relativePosition_relativePosition_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngles_frame_b_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder1_r_0_3_(fixed=false, unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngles_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Inertia bodyCylinder_body_I_32(fixed=false, min=-1E+060)
    " (3,2) element of inertia tensor";
    parameter Modelica_SIunits_Inertia bodyCylinder_body_I_31(fixed=false, min=-1E+060)
    " (3,1) element of inertia tensor";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativePosition_zeroPosition_frame_resolve_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Distance revolute2_cylinderLength(fixed=false)
    "Length of cylinder representing the joint axis";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativePosition_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Integer bodyCylinder1_body_cylinderColor_3_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Boolean rev_useAxisFlange(fixed=false)
    "= true, if axis flange is enabled";
    parameter Real bodyCylinder1_body_Q_start_3_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Modelica_SIunits_Distance bodyShape_height(fixed=false)
    " Height of shape.";
    parameter Real bodyCylinder_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyShape_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_v_0_2_(fixed=false, unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    parameter Modelica_SIunits_Diameter body_sphereDiameter(fixed=false)
    "Diameter of sphere";
    parameter Real bodyCylinder1_body_R_start_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real world_frame_b_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real world_frame_b_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_I_2_3_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativePosition_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyShape_lengthDirection_1_(fixed=false, unit="1")
    " Vector in length direction of shape, resolved in frame_a";
    parameter Real relativeAngles1_R_rel_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_frame_b_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Integer bodyShape_sequence_angleStates_1_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real revolute2_R_rel_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngles1_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real prismatic_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_w_a_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of frame_a resolved in frame_a";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer revolute2_cylinderColor_2_(max=255, fixed=false, min=0)
    "Color of cylinder representing the joint axis";
    parameter Real bodyShape_body_z_a_start_1_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_a = der(w_a), i.e., time derivative of angular velocity resolved in frame_a";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_R1_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_I_2_2_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real bodyShape_widthDirection_3_(fixed=false, unit="1")
    " Vector in width direction of shape, resolved in frame_a";
    parameter Real bodyCylinder1_frameTranslation_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder1_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real der_bodyShape_v_0_2_(fixed=false, unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_R1_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real prismatic_frame_b_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer world_gravitySphereColor_1_(max=255, fixed=false, min=0)
    "Color of gravity sphere";
    parameter Integer world_axisColor_y_2_(max=255, fixed=false, min=0);
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Modelica_SIunits_Inertia bodyCylinder_body_I_33(fixed=false, min=0.0)
    " (3,3) element of inertia tensor";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Boolean bodyCylinder_body_enforceStates(fixed=false)
    " = true, if absolute variables of body object shall be used as states (StateSelect.always)";
    parameter Real bodyCylinder1_frameTranslation_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyCylinder1_frameTranslation_color_3_(max=255, fixed=false, min=0)
    " Color of shape";
    parameter Real bodyShape_frameTranslation_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_b_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder1_body_I_3_2_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Modelica_SIunits_Length bodyCylinder_frameTranslation_length(fixed=false)
    " Length of shape";
    parameter Real body_Q_start_1_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Integer world_axisColor_x_2_(max=255, fixed=false, min=0)
    "Color of x-arrow";
    parameter Real bodyShape_body_w_0_start_2_(fixed=false, unit="rad/s")
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Integer bodyCylinder_body_sequence_angleStates_3_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Integer bodyShape_sequence_start_3_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Integer world_gravityArrowColor_3_(max=255, fixed=false, min=0)
    "Color of gravity arrow";
    parameter Boolean bodyCylinder_animation(fixed=false)
    "= true, if animation shall be enabled (show cylinder between frame_a and frame_b)";
    parameter Real relativePosition_zeroPosition_frame_resolve_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real revolute2_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_b_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_frameTranslation_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder1_body_w_0_start_2_(fixed=false, unit="rad/s")
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real bodyCylinder_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_frameTranslation_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_r_rel_3_(fixed=false, unit="m")
    "Relative position vector resolved in frame defined by resolveInFrame";
    parameter Real bodyCylinder_frameTranslation_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real prismatic_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_frameTranslation_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real rev_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_SIunits_Distance bodyCylinder1_frameTranslation_height(fixed=false)
    " Height of shape.";
    parameter Real relativeAngles_frame_a_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_I_3_1_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real relativeVelocity_relativePosition_relativePosition_r_rel_2_(fixed=false, unit="m")
    "Relative position vector frame_b.r_0 - frame_a.r_0 resolved in frame defined by resolveInFrame";
    parameter Real relativeAngularVelocity_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real der_bodyCylinder1_v_0_3_(fixed=false, unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    parameter Real rev_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean world_enableAnimation(fixed=false)
    "= true, if animation of all components is enabled";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_g_0_2_(fixed=false, unit="m/s2")
    "Gravity acceleration resolved in world frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyCylinder1_body_cylinderColor_2_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real bodyCylinder1_frameTranslation_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real body_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_body_g_0_1_(fixed=false, unit="m/s2")
    "Gravity acceleration resolved in world frame";
    parameter Real bodyCylinder1_frameTranslation_lengthDirection_1_(fixed=false, unit="1")
    " Vector in length direction of shape, resolved in frame_a";
    parameter Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB relativeAngularVelocity_relativeAngularVelocity_resolveInFrame(fixed=false)
    "Frame in which output vector w_rel is resolved (1: world, 2: frame_a, 3: frame_b, 4: frame_resolve)";
    parameter Integer bodyShape_body_sphereColor_2_(max=255, fixed=false, min=0)
    "Color of sphere";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Integer bodyShape_body_sequence_start_2_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Real bodyCylinder_body_R_start_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_r_0_3_(fixed=false, unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    parameter Boolean bodyShape_body_w_0_fixed(fixed=false)
    "= true, if w_0_start are used as initial values, else as guess values";
    parameter Real rev_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_R_start_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real rev_e_2_(fixed=false, unit="1")
    "Unit vector in direction of rotation axis, resolved in frame_a (= same as in frame_b)";
    parameter Real bodyCylinder_I_1_1_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean bodyCylinder_body_w_0_fixed(fixed=false)
    "= true, if w_0_start are used as initial values, else as guess values";
    parameter Boolean bodyCylinder1_body_enforceStates(fixed=false)
    " = true, if absolute variables of body object shall be used as states (StateSelect.always)";
    parameter Integer body_sequence_start_1_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Real bodyCylinder1_r_2_(fixed=false, unit="m")
    "Vector from frame_a to frame_b, resolved in frame_a";
    parameter Real bodyCylinder_body_v_0_2_(fixed=false, unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    parameter Real relativeVelocity_relativePosition_frame_b_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean bodyCylinder_w_0_fixed(fixed=false)
    "= true, if w_0_start are used as initial values, else as guess values";
    parameter Real bodyCylinder_frameTranslation_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Integer bodyCylinder1_body_sphereColor_1_(max=255, fixed=false, min=0)
    "Color of sphere";
    parameter Real body_Q_4_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Real der_relativeVelocity_relativePosition_relativePosition_frame_b_r_0_3_(fixed=false, unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    parameter Real bodyCylinder_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real prismatic_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean body_animation(fixed=false)
    "= true, if animation shall be enabled (show cylinder and sphere)";
    parameter Real bodyCylinder1_r_1_(fixed=false, unit="m")
    "Vector from frame_a to frame_b, resolved in frame_a";
    parameter Integer bodyCylinder1_body_sequence_start_2_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter StateSelect prismatic_stateSelect(fixed=false)
    "Priority to use distance s and v=der(s) as states";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_w_rel_2_(fixed=false, unit="rad/s")
    "Relative angular velocity vector";
    parameter Integer bodyCylinder_body_sequence_start_3_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Real bodyShape_r_1_(fixed=false, unit="m")
    "Vector from frame_a to frame_b resolved in frame_a";
    parameter Modelica_SIunits_Length world_lineWidth(fixed=false);
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real revolute2_R_rel_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer bodyShape_body_cylinderColor_1_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_phi_3_(displayUnit="deg", unit="rad", fixed=false)
    "Dummy or 3 angles to rotate world frame into frame_a of body";
    parameter Real relativePosition_relativePosition_frame_b_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativePosition_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder1_body_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder_body_Q_start_2_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Real body_I_3_3_(fixed=false, unit="kg.m2") "inertia tensor";
    parameter Real relativePosition_relativePosition_frame_resolve_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_R1_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real world_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_body_I_2_1_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Boolean bodyCylinder_angles_fixed(fixed=false)
    "= true, if angles_start are used as initial values, else as guess values";
    parameter Real relativeAngles1_R_rel_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_R1_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real prismatic_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_lengthDirection_2_(fixed=false, unit="1")
    " Vector in length direction of shape, resolved in frame_a";
    parameter Real rev_R_rel_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_SIunits_Length world_defaultForceLength(fixed=false)
    "Default for the fixed length of a shape representing a force (e.g., damper)";
    parameter Real relativeAngularVelocity1_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_SIunits_Distance bodyCylinder_radius(fixed=false)
    "Radius of cylinder";
    parameter Modelica_SIunits_Length bodyCylinder1_frameTranslation_length(fixed=false)
    " Length of shape";
    parameter Integer bodyCylinder1_color_3_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real relativeVelocity_frame_b_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngles_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_body_Q_1_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Real bodyCylinder1_body_angles_start_2_(displayUnit="deg", unit="rad", fixed=false)
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativePosition_relativePosition_frame_resolve_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Modelica_SIunits_Diameter bodyCylinder_body_sphereDiameter(fixed=false)
    "Diameter of sphere";
    parameter Real rev_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngles_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real rev_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real body_R_start_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_Mechanics_MultiBody_Types_SpecularCoefficient bodyCylinder_specularCoefficient(fixed=false)
    "Reflection of ambient light (= 0: light is completely absorbed)";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_frameTranslation_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real der_bodyShape_v_0_3_(fixed=false, unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    parameter Integer world_axisColor_y_3_(max=255, fixed=false, min=0);
    parameter Real bodyCylinder1_r_CM_3_(fixed=false, unit="m")
    "Position vector from frame_a to center of mass, resolved in frame_a";
    parameter Real body_g_0_2_(fixed=false, unit="m/s2")
    "Gravity acceleration resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_w_0_start_1_(fixed=false, unit="rad/s")
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frameTranslation_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyShape_body_sequence_start_3_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Real rev_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_body_phi_d_1_(fixed=false, unit="rad/s")
    "= der(phi)";
    parameter Real relativeAngles_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_w_rel_1_(fixed=false, unit="rad/s")
    "Relative angular velocity vector";
    parameter Real rev_e_3_(fixed=false, unit="1")
    "Unit vector in direction of rotation axis, resolved in frame_a (= same as in frame_b)";
    parameter Real bodyCylinder_r_shape_3_(fixed=false, unit="m")
    "Vector from frame_a to cylinder origin, resolved in frame_a";
    parameter Real bodyCylinder1_body_z_a_1_(fixed=false, unit="rad/s2")
    "Absolute angular acceleration of frame_a resolved in frame_a";
    parameter Integer bodyShape_sphereColor_3_(max=255, fixed=false, min=0)
    " Color of sphere of mass";
    parameter Real bodyCylinder1_body_phi_1_(displayUnit="deg", unit="rad", fixed=false)
    "Dummy or 3 angles to rotate world frame into frame_a of body";
    parameter Real body_I_1_3_(fixed=false, unit="kg.m2") "inertia tensor";
    parameter Real bodyCylinder_body_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Boolean bodyCylinder1_frameTranslation_animation(fixed=false)
    "= true, if animation shall be enabled";
    parameter Real relativePosition_relativePosition_frame_b_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_v_0_3_(fixed=false, unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    parameter Real bodyCylinder_body_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyCylinder1_color_2_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Integer bodyShape_sequence_start_2_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_body_z_a_2_(fixed=false, unit="rad/s2")
    "Absolute angular acceleration of frame_a resolved in frame_a";
    parameter Real relativePosition_zeroPosition_frame_resolve_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_frame_b_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder1_a_0_3_(fixed=false, unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    parameter Real relativePosition_relativePosition_frame_resolve_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Modelica_Mechanics_MultiBody_Types_ShapeExtra bodyCylinder_frameTranslation_extra(fixed=false, unit="1")
    " Additional parameter depending on shapeType (see docu of Visualizers.Advanced.Shape).";
    parameter Real relativeAngularVelocity_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real revolute2_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frameTranslation_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real der_relativeVelocity_der_r_rel_2__u(fixed=false)
    "der(Connector of Real input signal)";
    parameter Real bodyCylinder1_frameTranslation_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_r_CM_2_(fixed=false, unit="m")
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Real relativeVelocity_relativePosition_r_rel_2_(fixed=false, unit="m")
    "Relative position vector resolved in frame defined by resolveInFrame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_frameTranslation_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyShape_body_a_0_2_(fixed=false, unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    parameter Real der_bodyCylinder1_body_w_a_2_(fixed=false, unit="rad/s2")
    "der(Absolute angular velocity of frame_a resolved in frame_a)";
    parameter Real relativeVelocity_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Boolean bodyCylinder1_body_w_0_fixed(fixed=false)
    "= true, if w_0_start are used as initial values, else as guess values";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_z_a_2_(fixed=false, unit="rad/s2")
    "Absolute angular acceleration of frame_a resolved in frame_a";
    parameter Real relativeAngles_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Modelica_SIunits_Angle revolute2_phi_offset(fixed=false)
    "Relative angle offset (angle = phi_offset + phi)";
    parameter Integer bodyCylinder_body_sequence_start_2_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Modelica_SIunits_Length world_defaultJointLength(fixed=false)
    "Default for the fixed length of a shape representing a joint";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real rev_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_Q_3_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Real relativeAngularVelocity_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frameTranslation_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer body_sequence_start_2_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_frame_b_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter StateSelect revolute2_stateSelect(fixed=false)
    "Priority to use joint angle phi and w=der(phi) as states";
    parameter Real bodyCylinder_r_0_3_(fixed=false, unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    parameter Real der_relativeVelocity_relativePosition_relativePosition_frame_b_r_0_2_(fixed=false, unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    parameter Real body_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_frameTranslation_widthDirection_2_(fixed=false, unit="1")
    " Vector in width direction of shape, resolved in frame_a";
    parameter Real relativePosition_relativePosition_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Integer bodyShape_body_sphereColor_3_(max=255, fixed=false, min=0)
    "Color of sphere";
    parameter Real bodyCylinder1_body_Q_start_1_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Real bodyCylinder1_body_phi_start_1_(displayUnit="deg", unit="rad", fixed=false)
    "Potential angle states at initial time";
    parameter Real body_R_start_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean bodyCylinder1_animation(fixed=false)
    "= true, if animation shall be enabled (show cylinder between frame_a and frame_b)";
    parameter Real relativePosition_frame_b_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Modelica_Blocks_Interfaces_RealOutput relativeVelocity_der_r_rel_3__y(fixed=false)
    "Connector of Real output signal";
    parameter Real bodyShape_lengthDirection_3_(fixed=false, unit="1")
    " Vector in length direction of shape, resolved in frame_a";
    parameter Real rev_R_rel_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_Blocks_Interfaces_RealInput relativeVelocity_der_r_rel_3__u(fixed=false)
    "Connector of Real input signal";
    parameter Real relativeAngles_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder_frameTranslation_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngles1_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Integer bodyShape_color_2_(max=255, fixed=false, min=0)
    " Color of shape";
    parameter Real relativeVelocity_relativePosition_frame_b_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_pi(fixed=false);
    parameter Real bodyCylinder_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB relativeVelocity_tansformRelativeVector_frame_r_out(fixed=false)
    "Frame in which vector r_in shall be resolved and provided as r_out (1: world, 2: frame_a, 3: frame_b, 4: frame_resolve)";
    parameter Real bodyShape_body_z_a_start_3_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_a = der(w_a), i.e., time derivative of angular velocity resolved in frame_a";
    parameter Real prismatic_frame_b_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_R1_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_I_2_2_(fixed=false, unit="kg.m2") "inertia tensor";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Mass bodyCylinder1_m(fixed=false)
    "Mass of cylinder";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer world_gravitySphereColor_3_(max=255, fixed=false, min=0)
    "Color of gravity sphere";
    parameter Real bodyCylinder_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_g_0_3_(fixed=false, unit="m/s2")
    "Gravity acceleration resolved in world frame";
    parameter Real relativePosition_relativePosition_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyCylinder_body_sequence_start_1_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Real bodyShape_frameTranslation_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles1_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_phi_2_(displayUnit="deg", unit="rad", fixed=false)
    "Dummy or 3 angles to rotate world frame into frame_a of body";
    parameter Real body_R_start_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativePosition_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_frame_a_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real body_I_3_2_(fixed=false, unit="kg.m2") "inertia tensor";
    parameter Real body_Q_start_3_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Boolean world_animateWorld(fixed=false)
    "= true, if world coordinate system shall be visualized";
    parameter Real body_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_frameTranslation_r_shape_3_(fixed=false, unit="m")
    " Vector from frame_a to shape origin, resolved in frame_a";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder_body_R_start_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer world_gravityArrowColor_1_(max=255, fixed=false, min=0)
    "Color of gravity arrow";
    parameter Real revolute2_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_w_a_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of frame_a resolved in frame_a";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real prismatic_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Boolean bodyCylinder1_useQuaternions(fixed=false)
    " = true, if quaternions shall be used as potential states otherwise use 3 angles as potential states";
    parameter Real bodyShape_frameTranslation_lengthDirection_1_(fixed=false, unit="1")
    " Vector in length direction of shape, resolved in frame_a";
    parameter Real bodyShape_body_I_1_3_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real bodyCylinder_body_Q_2_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB relativeAngularVelocity1_resolveInFrame(fixed=false)
    "Frame in which output vector w_rel shall be resolved (1: world, 2: frame_a, 3: frame_b, 4: frame_resolve)";
    parameter Real rev_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder1_r_shape_3_(fixed=false, unit="m")
    "Vector from frame_a to cylinder origin, resolved in frame_a";
    parameter Real relativePosition_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder1_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngles_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real rev_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_frame_b_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_frameTranslation_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real rev_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_Q_start_4_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_Q_2_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Boolean bodyShape_frameTranslation_animation(fixed=false)
    "= true, if animation shall be enabled";
    parameter Integer body_sequence_start_3_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Integer bodyShape_body_sequence_angleStates_1_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real relativeVelocity_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_frame_b_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_r_0_2_(fixed=false, unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_body_phi_start_3_(displayUnit="deg", unit="rad", fixed=false)
    "Potential angle states at initial time";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_widthDirection_3_(fixed=false, unit="1")
    " Vector in width direction of shape, resolved in frame_a";
    parameter Real relativePosition_relativePosition_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real rev_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_w_rel_1_(fixed=false, unit="1/s")
    "Relative angular velocity vector between frame_a and frame_b resolved in frame defined by resolveInFrame";
    parameter Real bodyCylinder1_body_R_start_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_R1_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Integer bodyShape_sphereColor_2_(max=255, fixed=false, min=0)
    " Color of sphere of mass";
    parameter Real relativeVelocity_relativePosition_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_R_rel_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_R_rel_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer body_sequence_angleStates_3_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real body_I_1_2_(fixed=false, unit="kg.m2") "inertia tensor";
    parameter Real relativePosition_relativePosition_frame_b_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Mass bodyCylinder1_mo(fixed=false)
    "Mass of cylinder without hole";
    parameter Modelica_SIunits_Mass bodyCylinder1_mi(fixed=false)
    "Mass of hole of cylinder";
    parameter Real bodyCylinder_body_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Integer bodyShape_color_3_(max=255, fixed=false, min=0)
    " Color of shape";
    parameter Real bodyCylinder1_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_r_3_(fixed=false, unit="m")
    "Vector from frame_a to frame_b resolved in frame_a";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_r_in_3_(fixed=false)
    "Input vector resolved in frame defined by frame_r_in";
    parameter Real relativePosition_relativePosition_r_rel_3_(fixed=false, unit="m")
    "Relative position vector frame_b.r_0 - frame_a.r_0 resolved in frame defined by resolveInFrame";
    parameter Real bodyCylinder1_body_R_start_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngles1_R_rel_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Boolean body_useQuaternions(fixed=false)
    " = true, if quaternions shall be used as potential states otherwise use 3 angles as potential states";
    parameter Real bodyShape_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_Mechanics_MultiBody_Types_SpecularCoefficient bodyCylinder1_frameTranslation_specularCoefficient(fixed=false)
    "Reflection of ambient light (= 0: light is completely absorbed)";
    parameter Real revolute2_e_1_(fixed=false, unit="1")
    "Unit vector in direction of rotation axis, resolved in frame_a (= same as in frame_b)";
    parameter Real bodyShape_body_z_a_start_2_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_a = der(w_a), i.e., time derivative of angular velocity resolved in frame_a";
    parameter Real bodyCylinder1_frameTranslation_lengthDirection_3_(fixed=false, unit="1")
    " Vector in length direction of shape, resolved in frame_a";
    parameter Real bodyShape_v_0_2_(fixed=false, unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    parameter Real bodyShape_body_a_0_3_(fixed=false, unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer bodyCylinder1_body_cylinderColor_1_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real relativeVelocity_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngles_frame_b_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder1_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real revolute2_R_rel_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyCylinder1_color_1_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real relativePosition_frame_b_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyShape_body_cylinderColor_3_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real relativeAngles1_R_rel_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB relativeVelocity_tansformRelativeVector_frame_r_in(fixed=false)
    "Frame in which vector r_in is resolved (1: world, 2: frame_a, 3: frame_b, 4: frame_resolve)";
    parameter Real body_Q_start_4_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Real body_I_3_1_(fixed=false, unit="kg.m2") "inertia tensor";
    parameter Real bodyCylinder1_I_2_1_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Diameter body_cylinderDiameter(fixed=false)
    "Diameter of cylinder";
    parameter Real bodyCylinder_frameTranslation_lengthDirection_2_(fixed=false, unit="1")
    " Vector in length direction of shape, resolved in frame_a";
    parameter Modelica_SIunits_Inertia bodyCylinder1_body_I_33(fixed=false, min=0.0)
    " (3,3) element of inertia tensor";
    parameter Modelica_SIunits_Inertia bodyCylinder1_body_I_31(fixed=false, min=-1E+060)
    " (3,1) element of inertia tensor";
    parameter Real relativeVelocity_relativePosition_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity1_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder1_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frameTranslation_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles1_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real revolute2_R_rel_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer relativeAngles_sequence_1_(max=3, fixed=false, min=1)
    "Angles are returned to rotate frame_a around axes sequence[1], sequence[2] and finally sequence[3] into frame_b";
    parameter Real body_z_a_start_3_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_a = der(w_a), i.e., time derivative of angular velocity resolved in frame_a";
    parameter Real relativeAngles1_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real prismatic_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngles_R_rel_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyCylinder1_frameTranslation_color_2_(max=255, fixed=false, min=0)
    " Color of shape";
    parameter StateSelect rev_stateSelect(fixed=false)
    "Priority to use joint angle phi and w=der(phi) as states";
    parameter Real body_R_start_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngles_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real world_n_2_(fixed=false, unit="1")
    "Direction of gravity resolved in world frame (gravity = g*n/length(n))";
    parameter Real relativeAngles1_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_I_1_2_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Distance bodyShape_width(fixed=false)
    " Width of shape";
    parameter Real relativeAngularVelocity1_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_Q_start_3_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_r_shape_2_(fixed=false, unit="m")
    "Vector from frame_a to cylinder origin, resolved in frame_a";
    parameter Real relativeVelocity_relativePosition_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_b_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real revolute2_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real rev_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngles_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_I_3_2_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real body_I_2_3_(fixed=false, unit="kg.m2") "inertia tensor";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_body_angles_start_1_(displayUnit="deg", unit="rad", fixed=false)
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Real bodyCylinder1_frameTranslation_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder_body_phi_dd_3_(fixed=false, unit="rad/s2")
    "= der(phi_d)";
    parameter Integer body_cylinderColor_3_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean bodyCylinder1_body_angles_fixed(fixed=false)
    "= true, if angles_start are used as initial values, else as guess values";
    parameter Real bodyCylinder1_body_R_start_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_phi_d_3_(fixed=false, unit="rad/s")
    "= der(phi)";
    parameter Integer bodyShape_sphereColor_1_(max=255, fixed=false, min=0)
    " Color of sphere of mass";
    parameter Real bodyCylinder_body_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real rev_e_1_(fixed=false, unit="1")
    "Unit vector in direction of rotation axis, resolved in frame_a (= same as in frame_b)";
    parameter Modelica_SIunits_Length world_lineLength(fixed=false);
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Integer body_sequence_angleStates_2_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real bodyCylinder1_body_a_0_3_(fixed=false, unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_r_out_2_(fixed=false)
    "Input vector r_in resolved in frame defined by frame_r_out";
    parameter Real bodyCylinder_r_shape_1_(fixed=false, unit="m")
    "Vector from frame_a to cylinder origin, resolved in frame_a";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_frameTranslation_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real rev_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_r_shape_3_(fixed=false, unit="m")
    " Vector from frame_a to shape origin, resolved in frame_a";
    parameter Real bodyShape_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter StateSelect damper1_stateSelect(fixed=false)
    "Priority to use phi_rel and w_rel as states";
    parameter Real bodyShape_frameTranslation_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_phi_dd_1_(fixed=false, unit="rad/s2")
    "= der(phi_d)";
    parameter Integer bodyShape_body_sphereColor_1_(max=255, fixed=false, min=0)
    "Color of sphere";
    parameter Real prismatic_frame_b_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_z_a_2_(fixed=false, unit="rad/s2")
    "Absolute angular acceleration of frame_a resolved in frame_a";
    parameter Real bodyCylinder_frameTranslation_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_z_a_1_(fixed=false, unit="rad/s2")
    "Absolute angular acceleration of frame_a resolved in frame_a";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_R1_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngles1_R_rel_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frameTranslation_r_1_(fixed=false, unit="m")
    "Vector from frame_a to frame_b resolved in frame_a";
    parameter Real der_bodyCylinder_r_0_3_(fixed=false, unit="m/s")
    "der(Position vector from origin of world frame to origin of frame_a)";
    parameter Integer bodyCylinder1_body_sequence_angleStates_1_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real relativeAngles1_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_Q_1_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_r_in_2_(fixed=false)
    "Input vector resolved in frame defined by frame_r_in";
    parameter Real bodyCylinder1_frameTranslation_lengthDirection_2_(fixed=false, unit="1")
    " Vector in length direction of shape, resolved in frame_a";
    parameter Real relativePosition_relativePosition_frame_b_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_phi_3_(displayUnit="deg", unit="rad", fixed=false)
    "Dummy or 3 angles to rotate world frame into frame_a of body";
    parameter Real bodyShape_body_z_0_start_1_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real bodyCylinder_body_phi_start_2_(displayUnit="deg", unit="rad", fixed=false)
    "Potential angle states at initial time";
    parameter Real relativeAngularVelocity1_frame_b_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer bodyShape_body_sequence_start_1_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Real bodyShape_frameTranslation_frame_b_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real rev_n_1_(fixed=false, unit="1")
    "Axis of rotation resolved in frame_a (= same as in frame_b)";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_z_a_start_2_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_a = der(w_a), i.e., time derivative of angular velocity resolved in frame_a";
    parameter Real relativeAngularVelocity1_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_v_0_3_(fixed=false, unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    parameter Real body_R_start_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real world_n_3_(fixed=false, unit="1")
    "Direction of gravity resolved in world frame (gravity = g*n/length(n))";
    parameter Integer world_ndim_pointGravity(fixed=false);
    parameter Real bodyCylinder1_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_r_2_(fixed=false, unit="m")
    "Vector from frame_a to frame_b resolved in frame_a";
    parameter Real body_Q_1_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Real bodyCylinder1_body_R_start_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder1_frameTranslation_r_shape_2_(fixed=false, unit="m")
    " Vector from frame_a to shape origin, resolved in frame_a";
    parameter Real bodyCylinder_body_angles_start_2_(displayUnit="deg", unit="rad", fixed=false)
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Modelica_SIunits_Inertia bodyCylinder1_body_I_21(fixed=false, min=-1E+060)
    " (2,1) element of inertia tensor";
    parameter Modelica_SIunits_Inertia bodyCylinder1_body_I_22(fixed=false, min=0.0)
    " (2,2) element of inertia tensor";
    parameter Real body_R_start_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Length force_s_support(fixed=false)
    "Absolute position of support flange";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real revolute2_R_rel_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyCylinder1_frameTranslation_color_1_(max=255, fixed=false, min=0)
    " Color of shape";
    parameter Integer bodyCylinder_color_2_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder1_body_phi_dd_2_(fixed=false, unit="rad/s2")
    "= der(phi_d)";
    parameter Integer body_cylinderColor_2_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real relativePosition_zeroPosition_frame_resolve_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_body_phi_dd_2_(fixed=false, unit="rad/s2")
    "= der(phi_d)";
    parameter Real body_R_start_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativePosition_frame_b_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyShape_body_cylinderColor_2_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_a_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_lengthDirection_3_(fixed=false, unit="1")
    " Vector in length direction of shape, resolved in frame_a";
    parameter Real relativePosition_r_rel_3_(fixed=false, unit="m")
    "Relative position vector resolved in frame defined by resolveInFrame";
    parameter Real bodyShape_frameTranslation_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_relativePosition_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_body_phi_2_(displayUnit="deg", unit="rad", fixed=false)
    "Dummy or 3 angles to rotate world frame into frame_a of body";
    parameter Real bodyShape_body_v_0_3_(fixed=false, unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    parameter Modelica_SIunits_Angle revolute2_constantTorque_phi_support(fixed=false)
    "Absolute angle of support flange";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder1_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_b_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_lengthDirection_1_(fixed=false, unit="1")
    "Vector in length direction of cylinder, resolved in frame_a";
    parameter Real bodyShape_frameTranslation_lengthDirection_3_(fixed=false, unit="1")
    " Vector in length direction of shape, resolved in frame_a";
    parameter Real bodyCylinder1_frameTranslation_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_SIunits_Inertia bodyCylinder_I22(fixed=false)
    "Inertia with respect to axis through center of mass, perpendicular to cylinder axis";
    parameter Real prismatic_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder1_I_1_1_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real bodyShape_body_I_3_1_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real bodyShape_frameTranslation_frame_b_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real der_bodyCylinder_v_0_3_(fixed=false, unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    parameter Real bodyCylinder1_body_z_0_start_1_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real bodyCylinder1_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_r_3_(fixed=false, unit="m")
    "Vector from frame_a to frame_b resolved in frame_a";
    parameter Modelica_SIunits_Distance bodyCylinder1_innerRadius(fixed=false)
    "Inner-Radius of cylinder";
    parameter Real relativeAngles_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_body_I_1_1_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real bodyShape_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real body_g_0_1_(fixed=false, unit="m/s2")
    "Gravity acceleration resolved in world frame";
    parameter Real bodyCylinder1_r_shape_1_(fixed=false, unit="m")
    "Vector from frame_a to cylinder origin, resolved in frame_a";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer prismatic_boxColor_2_(max=255, fixed=false, min=0)
    "Color of prismatic joint box";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_r_out_3_(fixed=false)
    "Input vector r_in resolved in frame defined by frame_r_out";
    parameter Real relativeAngles1_frame_b_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real revolute2_n_2_(fixed=false, unit="1")
    "Axis of rotation resolved in frame_a (= same as in frame_b)";
    parameter Real relativePosition_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_b_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity1_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Integer bodyShape_body_sequence_angleStates_3_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real bodyCylinder1_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real revolute2_e_3_(fixed=false, unit="1")
    "Unit vector in direction of rotation axis, resolved in frame_a (= same as in frame_b)";
    parameter Real bodyCylinder_frameTranslation_widthDirection_1_(fixed=false, unit="1")
    " Vector in width direction of shape, resolved in frame_a";
    parameter Real rev_R_rel_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean body_w_0_fixed(fixed=false)
    "= true, if w_0_start are used as initial values, else as guess values";
    parameter Real bodyCylinder1_frameTranslation_widthDirection_2_(fixed=false, unit="1")
    " Vector in width direction of shape, resolved in frame_a";
    parameter Real bodyCylinder_body_z_a_start_3_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_a = der(w_a), i.e., time derivative of angular velocity resolved in frame_a";
    parameter Real rev_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_body_r_0_3_(fixed=false, unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    parameter Real bodyShape_body_z_0_start_2_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real body_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_frameTranslation_frame_b_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_body_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_body_Q_start_2_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Integer bodyCylinder_body_sphereColor_1_(max=255, fixed=false, min=0)
    "Color of sphere";
    parameter Real der_bodyCylinder1_r_0_3_(fixed=false, unit="m/s")
    "der(Position vector from origin of world frame to origin of frame_a)";
    parameter Real bodyCylinder_body_R_start_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_phi_d_2_(fixed=false, unit="rad/s")
    "= der(phi)";
    parameter Modelica_Mechanics_MultiBody_Types_ShapeType bodyCylinder_frameTranslation_shapeType =  " Type of shape";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_R1_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_SIunits_Length world_defaultArrowDiameter(fixed=false)
    "Default for arrow diameter (e.g., of forces, torques, sensors)";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativePosition_frame_a_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_body_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer body_sequence_angleStates_1_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real revolute2_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean bodyShape_body_enforceStates(fixed=false)
    " = true, if absolute variables of body object shall be used as states (StateSelect.always)";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyShape_color_1_(max=255, fixed=false, min=0)
    " Color of shape";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_body_phi_start_2_(displayUnit="deg", unit="rad", fixed=false)
    "Potential angle states at initial time";
    parameter Real bodyCylinder_body_z_0_start_3_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real relativePosition_relativePosition_frame_resolve_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_frame_b_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_body_I_3_3_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real prismatic_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frameTranslation_r_shape_2_(fixed=false, unit="m")
    " Vector from frame_a to shape origin, resolved in frame_a";
    parameter Integer bodyShape_frameTranslation_color_3_(max=255, fixed=false, min=0)
    " Color of shape";
    parameter Integer relativeAngles1_sequence_1_(max=3, fixed=false, min=1)
    "Angles are returned to rotate frame_a around axes sequence[1], sequence[2] and finally sequence[3] into frame_b";
    parameter Real der_bodyCylinder_r_0_2_(fixed=false, unit="m/s")
    "der(Position vector from origin of world frame to origin of frame_a)";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real body_z_a_1_(fixed=false, unit="rad/s2")
    "Absolute angular acceleration of frame_a resolved in frame_a";
    parameter Real bodyShape_frameTranslation_r_2_(fixed=false, unit="m")
    "Vector from frame_a to frame_b resolved in frame_a";
    parameter Real relativeAngularVelocity1_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB relativePosition_relativePosition_resolveInFrame(fixed=false)
    "Frame in which output vector r_rel is resolved (1: world, 2: frame_a, 3: frame_b, 4: frame_resolve)";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeAngles1_frame_a_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real body_R_start_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real prismatic_n_2_(fixed=false, unit="1")
    "Axis of translation resolved in frame_a (= same as in frame_b)";
    parameter Real bodyShape_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder1_I_2_3_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_body_g_0_1_(fixed=false, unit="m/s2")
    "Gravity acceleration resolved in world frame";
    parameter Modelica_SIunits_Inertia bodyCylinder1_body_I_11(fixed=false, min=0.0)
    " (1,1) element of inertia tensor";
    parameter Real relativePosition_r_rel_2_(fixed=false, unit="m")
    "Relative position vector resolved in frame defined by resolveInFrame";
    parameter Real relativeAngles_R_rel_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frame_b_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Length world_defaultJointWidth(fixed=false)
    "Default for the fixed width of a shape representing a joint";
    parameter Real relativeAngles1_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_v_rel_3_(fixed=false, unit="m/s")
    "Relative velocity vector resolved in frame defined by resolveInFrame";
    parameter Real body_z_a_start_1_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_a = der(w_a), i.e., time derivative of angular velocity resolved in frame_a";
    parameter Real bodyShape_frameTranslation_frame_b_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_frame_b_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer relativeAngles_sequence_3_(max=3, fixed=false, min=1)
    "Angles are returned to rotate frame_a around axes sequence[1], sequence[2] and finally sequence[3] into frame_b";
    parameter Real bodyShape_body_phi_d_2_(fixed=false, unit="rad/s")
    "= der(phi)";
    parameter Real relativeAngles_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngles1_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer body_cylinderColor_1_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real bodyCylinder_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_frame_b_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_r_in_3_(fixed=false)
    "Input vector resolved in frame defined by frame_r_in";
    parameter Real relativeVelocity_relativePosition_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean damper_useHeatPort(fixed=false)
    "=true, if heatPort is enabled";
    parameter Integer prismatic_boxColor_3_(max=255, fixed=false, min=0)
    "Color of prismatic joint box";
    parameter Real relativeAngles1_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_frame_b_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngles_R_rel_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Angle rev_phi_offset(fixed=false)
    "Relative angle offset (angle = phi_offset + phi)";
    parameter Real bodyCylinder_body_Q_3_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Modelica_SIunits_Mass bodyCylinder_mi(fixed=false)
    "Mass of hole of cylinder";
    parameter Real bodyCylinder1_I_3_1_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real bodyCylinder_body_I_1_3_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real bodyShape_body_angles_start_3_(displayUnit="deg", unit="rad", fixed=false)
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Real relativeAngles_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_phi_start_1_(displayUnit="deg", unit="rad", fixed=false)
    "Potential angle states at initial time";
    parameter Real body_phi_start_1_(displayUnit="deg", unit="rad", fixed=false)
    "Potential angle states at initial time";
    parameter Real bodyCylinder_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Boolean revolute2_useAxisFlange(fixed=false)
    "= true, if axis flange is enabled";
    parameter Real body_I_2_1_(fixed=false, unit="kg.m2") "inertia tensor";
    parameter Real relativePosition_zeroPosition_frame_resolve_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_lengthDirection_3_(fixed=false, unit="1")
    "Vector in length direction of cylinder, resolved in frame_a";
    parameter Real bodyCylinder1_body_R_start_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_Q_start_3_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Boolean body_z_0_fixed(fixed=false)
    "= true, if z_0_start are used as initial values, else as guess values";
    parameter Real relativeVelocity_relativePosition_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder1_body_I_1_1_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Modelica_SIunits_Inertia bodyShape_body_I_33(fixed=false, min=0.0)
    " (3,3) element of inertia tensor";
    parameter Modelica_SIunits_Inertia bodyShape_body_I_32(fixed=false, min=-1E+060)
    " (3,2) element of inertia tensor";
    parameter Modelica_SIunits_Inertia bodyShape_body_I_31(fixed=false, min=-1E+060)
    " (3,1) element of inertia tensor";
    parameter Real bodyShape_body_v_0_2_(fixed=false, unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    parameter Real bodyShape_frameTranslation_widthDirection_3_(fixed=false, unit="1")
    " Vector in width direction of shape, resolved in frame_a";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_r_0_2_(fixed=false, unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frameTranslation_lengthDirection_2_(fixed=false, unit="1")
    " Vector in length direction of shape, resolved in frame_a";
    parameter Real bodyCylinder1_frameTranslation_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real der_bodyCylinder_v_0_2_(fixed=false, unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    parameter Real bodyShape_frameTranslation_frame_b_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_body_r_0_2_(fixed=false, unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    parameter Real bodyCylinder1_body_Q_3_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Real bodyShape_frameTranslation_r_shape_3_(fixed=false, unit="m")
    " Vector from frame_a to shape origin, resolved in frame_a";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyShape_frameTranslation_color_2_(max=255, fixed=false, min=0)
    " Color of shape";
    parameter Real relativeVelocity_relativePosition_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_R_start_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngles1_R_rel_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Integer bodyCylinder1_body_sequence_angleStates_3_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Integer bodyShape_body_sequence_angleStates_2_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real bodyCylinder1_body_phi_dd_3_(fixed=false, unit="rad/s2")
    "= der(phi_d)";
    parameter Real relativeAngularVelocity1_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder1_frameTranslation_widthDirection_1_(fixed=false, unit="1")
    " Vector in width direction of shape, resolved in frame_a";
    parameter Integer bodyCylinder_frameTranslation_color_1_(max=255, fixed=false, min=0)
    " Color of shape";
    parameter Real bodyShape_frameTranslation_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real revolute2_e_2_(fixed=false, unit="1")
    "Unit vector in direction of rotation axis, resolved in frame_a (= same as in frame_b)";
    parameter Real bodyShape_body_z_0_start_3_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real rev_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real rev_R_rel_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_r_shape_2_(fixed=false, unit="m")
    "Vector from frame_a to cylinder origin, resolved in frame_a";
    parameter Real relativeAngles_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Boolean bodyShape_useQuaternions(fixed=false)
    " = true, if quaternions shall be used as potential states otherwise use 3 angles as potential states";
    parameter Real bodyShape_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_v_rel_2_(fixed=false, unit="m/s")
    "Relative velocity vector resolved in frame defined by resolveInFrame";
    parameter Real bodyCylinder1_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frameTranslation_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real der_bodyCylinder_body_v_0_2_(fixed=false, unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    parameter Real relativeAngularVelocity_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real world_n_1_(fixed=false, unit="1")
    "Direction of gravity resolved in world frame (gravity = g*n/length(n))";
    parameter Real bodyCylinder_body_phi_1_(displayUnit="deg", unit="rad", fixed=false)
    "Dummy or 3 angles to rotate world frame into frame_a of body";
    parameter Real revolute2_R_rel_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_g_0_2_(fixed=false, unit="m/s2")
    "Gravity acceleration resolved in world frame";
    parameter Real bodyShape_body_phi_d_3_(fixed=false, unit="rad/s")
    "= der(phi)";
    parameter Real rev_n_3_(fixed=false, unit="1")
    "Axis of rotation resolved in frame_a (= same as in frame_b)";
    parameter Real bodyCylinder_frameTranslation_r_2_(fixed=false, unit="m")
    "Vector from frame_a to frame_b resolved in frame_a";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real revolute2_n_3_(fixed=false, unit="1")
    "Axis of rotation resolved in frame_a (= same as in frame_b)";
    parameter Modelica_SIunits_Inertia bodyCylinder1_body_I_32(fixed=false, min=-1E+060)
    " (3,2) element of inertia tensor";
    parameter Real body_Q_3_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Real bodyCylinder_frameTranslation_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_z_0_start_2_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Boolean bodyShape_body_animation(fixed=false)
    "= true, if animation shall be enabled (show cylinder and sphere)";
    parameter Real body_R_start_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_Blocks_Interfaces_RealInput relativeVelocity_der_r_rel_2__u(fixed=false)
    "Connector of Real input signal";
    parameter Modelica_Blocks_Interfaces_RealOutput relativeVelocity_der_r_rel_2__y(fixed=false)
    "Connector of Real output signal";
    parameter Real relativeAngles_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_z_a_start_3_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_a = der(w_a), i.e., time derivative of angular velocity resolved in frame_a";
    parameter Real bodyShape_body_angles_start_2_(displayUnit="deg", unit="rad", fixed=false)
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Real bodyCylinder1_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder1_I_3_2_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real bodyCylinder1_body_v_0_3_(fixed=false, unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    parameter Real body_phi_start_2_(displayUnit="deg", unit="rad", fixed=false)
    "Potential angle states at initial time";
    parameter Real bodyShape_body_z_a_1_(fixed=false, unit="rad/s2")
    "Absolute angular acceleration of frame_a resolved in frame_a";
    parameter Real world_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_lengthDirection_1_(fixed=false, unit="1")
    " Vector in length direction of shape, resolved in frame_a";
    parameter Real bodyCylinder_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real prismatic_n_3_(fixed=false, unit="1")
    "Axis of translation resolved in frame_a (= same as in frame_b)";
    parameter Real relativePosition_relativePosition_frame_resolve_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyShape_body_g_0_2_(fixed=false, unit="m/s2")
    "Gravity acceleration resolved in world frame";
    parameter Real relativeVelocity_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_body_r_CM_1_(fixed=false, unit="m")
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Real bodyCylinder1_I_2_2_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real relativeAngularVelocity_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity1_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder_a_0_2_(fixed=false, unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    parameter Boolean revolute2_animation(fixed=false)
    "= true, if animation shall be enabled (show axis as cylinder)";
    parameter Real bodyCylinder1_frameTranslation_r_2_(fixed=false, unit="m")
    "Vector from frame_a to frame_b resolved in frame_a";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder1_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frameTranslation_widthDirection_2_(fixed=false, unit="1")
    " Vector in width direction of shape, resolved in frame_a";
    parameter Real bodyShape_body_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB relativeVelocity_resolveInFrame(fixed=false)
    "Frame in which output vector v_rel shall be resolved (1: world, 2: frame_a, 3: frame_b, 4: frame_resolve)";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_body_r_0_3_(fixed=false, unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_Mechanics_MultiBody_Types_SpecularCoefficient bodyCylinder1_specularCoefficient(fixed=false)
    "Reflection of ambient light (= 0: light is completely absorbed)";
    parameter Modelica_SIunits_Distance bodyShape_frameTranslation_width(fixed=false)
    " Width of shape";
    parameter Real bodyShape_body_I_3_3_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real bodyCylinder_lengthDirection_3_(fixed=false, unit="1")
    "Vector in length direction of cylinder, resolved in frame_a";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_v_0_3_(fixed=false, unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    parameter Real bodyCylinder_body_R_start_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Length world_nominalLength(fixed=false);
    parameter Real relativeAngles1_frame_b_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_r_in_2_(fixed=false)
    "Input vector resolved in frame defined by frame_r_in";
    parameter Real bodyShape_body_phi_start_3_(displayUnit="deg", unit="rad", fixed=false)
    "Potential angle states at initial time";
    parameter Real prismatic_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativePosition_frame_b_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativePosition_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Modelica_SIunits_Torque revolute2_fixed_flange_tau(fixed=false)
    "Cut torque in the flange";
    parameter Real relativeVelocity_relativePosition_frame_b_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_body_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Inertia bodyShape_body_I_21(fixed=false, min=-1E+060)
    " (2,1) element of inertia tensor";
    parameter Modelica_SIunits_Inertia bodyShape_body_I_22(fixed=false, min=0.0)
    " (2,2) element of inertia tensor";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_g_0_1_(fixed=false, unit="m/s2")
    "Gravity acceleration resolved in world frame";
    parameter Real relativeAngles_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real rev_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer bodyCylinder_frameTranslation_color_2_(max=255, fixed=false, min=0)
    " Color of shape";
    parameter Real bodyCylinder1_frameTranslation_r_1_(fixed=false, unit="m")
    "Vector from frame_a to frame_b resolved in frame_a";
    parameter Real rev_R_rel_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frame_b_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_r_shape_2_(fixed=false, unit="m")
    " Vector from frame_a to shape origin, resolved in frame_a";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_body_I_2_1_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real bodyCylinder_body_R_start_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Inertia bodyCylinder1_I22(fixed=false)
    "Inertia with respect to axis through center of mass, perpendicular to cylinder axis";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_I_1_3_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real revolute2_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_w_rel_1_(fixed=false, unit="1/s")
    "Relative angular velocity vector between frame_a and frame_b resolved in frame defined by resolveInFrame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_body_angles_start_1_(displayUnit="deg", unit="rad", fixed=false)
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_w_rel_2_(fixed=false, unit="rad/s")
    "Relative angular velocity vector";
    parameter Real relativePosition_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_Q_2_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Real der_bodyCylinder1_body_frame_a_r_0_3_(fixed=false, unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    parameter Real der_bodyShape_r_0_2_(fixed=false, unit="m/s")
    "der(Position vector from origin of world frame to origin of frame_a)";
    parameter Real bodyCylinder1_frameTranslation_r_shape_1_(fixed=false, unit="m")
    " Vector from frame_a to shape origin, resolved in frame_a";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_R_start_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_z_0_start_1_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real bodyShape_body_R_start_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real rev_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_I_3_1_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real relativePosition_relativePosition_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder1_body_Q_2_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Real prismatic_boxWidthDirection_1_(fixed=false, unit="1")
    "Vector in width direction of box, resolved in frame_a";
    parameter Integer bodyShape_frameTranslation_color_1_(max=255, fixed=false, min=0)
    " Color of shape";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_Mechanics_MultiBody_Types_SpecularCoefficient bodyCylinder_frameTranslation_specularCoefficient(fixed=false)
    "Reflection of ambient light (= 0: light is completely absorbed)";
    parameter Integer bodyCylinder1_sequence_start_1_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Real bodyCylinder1_lengthDirection_2_(fixed=false, unit="1")
    "Vector in length direction of cylinder, resolved in frame_a";
    parameter Real bodyShape_body_R_start_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Integer bodyCylinder1_body_sequence_angleStates_2_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Integer revolute2_cylinderColor_1_(max=255, fixed=false, min=0)
    "Color of cylinder representing the joint axis";
    parameter Real relativeAngles1_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeAngularVelocity_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_relativePosition_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_body_R_start_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_body_g_0_3_(fixed=false, unit="m/s2")
    "Gravity acceleration resolved in world frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_frameTranslation_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_body_R_start_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_r_CM_1_(fixed=false, unit="m")
    "Position vector from frame_a to center of mass, resolved in frame_a";
    parameter Real prismatic_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frameTranslation_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real der_bodyCylinder_body_v_0_3_(fixed=false, unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    parameter Boolean bodyShape_animation(fixed=false)
    "= true, if animation shall be enabled (show shape between frame_a and frame_b and optionally a sphere at the center of mass)";
    parameter Real rev_n_2_(fixed=false, unit="1")
    "Axis of rotation resolved in frame_a (= same as in frame_b)";
    parameter Real bodyCylinder_frameTranslation_r_1_(fixed=false, unit="m")
    "Vector from frame_a to frame_b resolved in frame_a";
    parameter Real der_body_v_0_3_(fixed=false, unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    parameter Real relativeVelocity_relativePosition_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Modelica_Mechanics_MultiBody_Types_SpecularCoefficient body_specularCoefficient(fixed=false)
    "Reflection of ambient light (= 0: light is completely absorbed)";
    parameter Real relativeAngles1_frame_b_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativePosition_frame_b_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativePosition_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles1_frame_b_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer prismatic_boxColor_1_(max=255, fixed=false, min=0)
    "Color of prismatic joint box";
    parameter Real der_bodyShape_body_frame_a_r_0_3_(fixed=false, unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    parameter Modelica_SIunits_Distance bodyCylinder_innerRadius(fixed=false)
    "Inner-Radius of cylinder";
    parameter Real relativeAngles_R_rel_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_z_a_start_2_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_a = der(w_a), i.e., time derivative of angular velocity resolved in frame_a";
    parameter Real relativeAngularVelocity1_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativePosition_zeroPosition_frame_resolve_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativePosition_relativePosition_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_I_3_3_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real bodyShape_body_z_a_2_(fixed=false, unit="rad/s2")
    "Absolute angular acceleration of frame_a resolved in frame_a";
    parameter Real bodyCylinder1_body_phi_start_3_(displayUnit="deg", unit="rad", fixed=false)
    "Potential angle states at initial time";
    parameter Real relativeAngles1_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_phi_start_3_(displayUnit="deg", unit="rad", fixed=false)
    "Potential angle states at initial time";
    parameter Integer bodyCylinder_body_sphereColor_2_(max=255, fixed=false, min=0)
    "Color of sphere";
    parameter Real bodyCylinder1_frameTranslation_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_body_R_start_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_zeroPosition_frame_resolve_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder1_body_w_a_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of frame_a resolved in frame_a";
    parameter Real bodyShape_body_Q_start_1_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Real bodyCylinder1_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyShape_body_I_2_2_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Boolean bodyCylinder_body_animation(fixed=false)
    "= true, if animation shall be enabled (show cylinder and sphere)";
    parameter Modelica_SIunits_Distance damper1_s_nominal(fixed=false)
    "Nominal value of s_rel (used for scaling)";
    parameter Real relativeAngularVelocity_w_rel_2_(fixed=false, unit="1/s")
    "Relative angular velocity vector between frame_a and frame_b resolved in frame defined by resolveInFrame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_relativePosition_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_phi_d_2_(fixed=false, unit="rad/s") "= der(phi)";
    parameter Integer world_ndim(fixed=false);
    parameter Real relativePosition_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder1_frameTranslation_r_3_(fixed=false, unit="m")
    "Vector from frame_a to frame_b resolved in frame_a";
    parameter Real relativeAngles1_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Inertia bodyShape_body_I_11(fixed=false, min=0.0)
    " (1,1) element of inertia tensor";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyShape_body_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_frameTranslation_widthDirection_1_(fixed=false, unit="1")
    " Vector in width direction of shape, resolved in frame_a";
    parameter Boolean bodyCylinder_frameTranslation_animation(fixed=false)
    "= true, if animation shall be enabled";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyCylinder_frameTranslation_color_3_(max=255, fixed=false, min=0)
    " Color of shape";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_lengthDirection_2_(fixed=false, unit="1")
    "Vector in length direction of cylinder, resolved in frame_a";
    parameter Real bodyShape_body_I_3_2_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real bodyCylinder_body_w_0_start_2_(fixed=false, unit="rad/s")
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_b_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Modelica_SIunits_Length world_scaledLabel(fixed=false);
    parameter Real bodyShape_frameTranslation_r_shape_1_(fixed=false, unit="m")
    " Vector from frame_a to shape origin, resolved in frame_a";
    parameter Real relativePosition_relativePosition_frame_a_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real prismatic_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_Q_1_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Real prismatic_boxWidthDirection_2_(fixed=false, unit="1")
    "Vector in width direction of box, resolved in frame_a";
    parameter Boolean bodyShape_w_0_fixed(fixed=false)
    "= true, if w_0_start are used as initial values, else as guess values";
    parameter Real bodyCylinder1_lengthDirection_1_(fixed=false, unit="1")
    "Vector in length direction of cylinder, resolved in frame_a";
    parameter Real revolute2_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_SIunits_Mass bodyCylinder_mo(fixed=false)
    "Mass of cylinder without hole";
    parameter Real bodyShape_frameTranslation_r_3_(fixed=false, unit="m")
    "Vector from frame_a to frame_b resolved in frame_a";
    parameter Real bodyShape_body_R_start_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder1_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Integer bodyCylinder1_body_sequence_start_3_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Real bodyCylinder_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_relativePosition_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Integer body_sphereColor_1_(max=255, fixed=false, min=0)
    "Color of sphere";
    parameter Real relativeAngles_frame_b_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder1_frameTranslation_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean bodyShape_animateSphere(fixed=false)
    "= true, if mass shall be animated as sphere provided animation=true";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB relativeAngularVelocity_resolveInFrame(fixed=false)
    "Frame in which output vector w_rel shall be resolved (1: world, 2: frame_a, 3: frame_b, 4: frame_resolve)";
    parameter Real bodyCylinder_frameTranslation_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Integer rev_cylinderColor_1_(max=255, fixed=false, min=0)
    "Color of cylinder representing the joint axis";
    parameter Modelica_SIunits_Angle phi1_start(fixed=false);
    parameter Real prismatic_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frame_b_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real world_frame_b_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_zeroPosition_frame_resolve_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_frameTranslation_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real revolute2_n_1_(fixed=false, unit="1")
    "Axis of rotation resolved in frame_a (= same as in frame_b)";
    parameter Real bodyShape_body_phi_d_1_(fixed=false, unit="rad/s")
    "= der(phi)";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles1_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_body_a_0_3_(fixed=false, unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real der_bodyShape_r_0_3_(fixed=false, unit="m/s")
    "der(Position vector from origin of world frame to origin of frame_a)";
    parameter Real body_Q_start_2_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB relativeVelocity_relativePosition_resolveInFrame(fixed=false)
    "Frame in which output vector r_rel shall be resolved (1: world, 2: frame_a, 3: frame_b, 4: frame_resolve)";
    parameter Real relativeAngularVelocity1_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Distance bodyCylinder1_frameTranslation_width(fixed=false)
    " Width of shape";
    parameter Real bodyShape_body_R_start_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real der_bodyShape_body_frame_a_r_0_2_(fixed=false, unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    parameter Real bodyShape_body_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles1_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_z_a_3_(fixed=false, unit="rad/s2")
    "Absolute angular acceleration of frame_a resolved in frame_a";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_frame_a_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder1_body_phi_d_2_(fixed=false, unit="rad/s")
    "= der(phi)";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real prismatic_n_1_(fixed=false, unit="1")
    "Axis of translation resolved in frame_a (= same as in frame_b)";
    parameter Modelica_SIunits_Length world_defaultAxisLength(fixed=false)
    "Default for length of a frame axis (but not world frame)";
    parameter Real bodyShape_frameTranslation_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_body_r_CM_3_(fixed=false, unit="m")
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Boolean force_useSupport(fixed=false)
    "= true, if support flange enabled, otherwise implicitly grounded";
    parameter Real bodyCylinder_body_R_start_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_SIunits_Mass body_m(fixed=false) "Mass of rigid body";
    parameter Real bodyShape_body_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_body_phi_3_(displayUnit="deg", unit="rad", fixed=false)
    "Dummy or 3 angles to rotate world frame into frame_a of body";
    parameter Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB relativePosition_resolveInFrame(fixed=false)
    "Frame in which output vector r_rel shall be resolved (1: world, 2: frame_a, 3: frame_b, 4: frame_resolve)";
    parameter Real bodyCylinder_body_R_start_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_R_start_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_pi(fixed=false);
    parameter Real relativeAngularVelocity1_frame_a_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_body_R_start_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Modelica_SIunits_Length bodyCylinder_length(fixed=false)
    "Length of cylinder";
    parameter Modelica_SIunits_Distance world_axisDiameter(fixed=false)
    "Diameter of world axes arrows";
    parameter Real relativeVelocity_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles1_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Length bodyShape_frameTranslation_length(fixed=false)
    " Length of shape";
    parameter Real bodyShape_body_phi_start_1_(displayUnit="deg", unit="rad", fixed=false)
    "Potential angle states at initial time";
    parameter Real relativePosition_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real prismatic_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_body_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyCylinder1_sequence_angleStates_2_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyShape_body_Q_4_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Real relativeAngles_frame_b_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativePosition_relativePosition_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_g_0_3_(fixed=false, unit="m/s2")
    "Gravity acceleration resolved in world frame";
    parameter Modelica_SIunits_Diameter bodyShape_body_cylinderDiameter(fixed=false)
    "Diameter of cylinder";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Modelica_SIunits_Distance bodyShape_frameTranslation_height(fixed=false)
    " Height of shape.";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean prismatic_useAxisFlange(fixed=false)
    "= true, if axis flange is enabled";
    parameter Modelica_SIunits_Mass bodyShape_m(fixed=false)
    "Mass of rigid body";
    parameter Integer body_sphereColor_2_(max=255, fixed=false, min=0)
    "Color of sphere";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Distance bodyCylinder_frameTranslation_height(fixed=false)
    " Height of shape.";
    parameter Real bodyShape_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_zeroPosition_frame_resolve_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real revolute2_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real prismatic_e_1_(fixed=false, unit="1")
    "Unit vector in direction of prismatic axis n";
    parameter Real bodyShape_body_I_2_3_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_body_angles_start_3_(displayUnit="deg", unit="rad", fixed=false)
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter StateSelect damper_stateSelect(fixed=false)
    "Priority to use phi_rel and w_rel as states";
    parameter Real body_phi_d_3_(fixed=false, unit="rad/s") "= der(phi)";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativePosition_relativePosition_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_resolve_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Modelica_Mechanics_MultiBody_Types_AxisLabel world_label1 =  "Label of horizontal axis in icon";
    parameter Modelica_Mechanics_MultiBody_Types_AxisLabel world_label2 =  "Label of vertical axis in icon";
    parameter Real relativePosition_relativePosition_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder1_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Modelica_SIunits_Length world_headLength(fixed=false);
    parameter Real rev_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real der_relativeVelocity_der_r_rel_3__u(fixed=false)
    "der(Connector of Real input signal)";
    parameter Integer relativeAngles_sequence_2_(max=3, fixed=false, min=1)
    "Angles are returned to rotate frame_a around axes sequence[1], sequence[2] and finally sequence[3] into frame_b";
    parameter Real relativePosition_relativePosition_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_b_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real rev_R_rel_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Modelica_SIunits_Diameter bodyCylinder1_body_cylinderDiameter(fixed=false)
    "Diameter of cylinder";
    parameter Real prismatic_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_b_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real prismatic_boxWidthDirection_3_(fixed=false, unit="1")
    "Vector in width direction of box, resolved in frame_a";
    parameter Real prismatic_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean bodyCylinder1_z_0_fixed(fixed=false)
    "= true, if z_0_start are used as initial values, else as guess values";
    parameter Integer bodyCylinder1_sequence_start_3_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Boolean bodyCylinder_body_useQuaternions(fixed=false)
    " = true, if quaternions shall be used as potential states otherwise use 3 angles as potential states";
    parameter Real revolute2_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity1_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_phi_d_3_(fixed=false, unit="rad/s")
    "= der(phi)";
    parameter Real bodyShape_frameTranslation_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_SIunits_Distance rev_cylinderLength(fixed=false)
    "Length of cylinder representing the joint axis";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_R_rel_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real pi(fixed=false);
    parameter Real bodyShape_body_r_CM_2_(fixed=false, unit="m")
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Real bodyCylinder_r_CM_3_(fixed=false, unit="m")
    "Position vector from frame_a to center of mass, resolved in frame_a";
    parameter Real bodyCylinder1_body_I_1_3_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real bodyShape_body_R_start_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer rev_cylinderColor_2_(max=255, fixed=false, min=0)
    "Color of cylinder representing the joint axis";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real prismatic_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_zeroPosition_frame_resolve_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_I_3_3_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Boolean bodyCylinder_z_0_fixed(fixed=false)
    "= true, if z_0_start are used as initial values, else as guess values";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_R_rel_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_z_a_start_1_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_a = der(w_a), i.e., time derivative of angular velocity resolved in frame_a";
    parameter Real bodyCylinder1_body_r_CM_1_(fixed=false, unit="m")
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Real relativePosition_relativePosition_frame_a_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder_body_Q_start_4_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_R_rel_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles1_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real revolute2_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Distance prismatic_boxHeight(fixed=false)
    "Height of prismatic joint box";
    parameter Boolean bodyCylinder1_body_useQuaternions(fixed=false)
    " = true, if quaternions shall be used as potential states otherwise use 3 angles as potential states";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Boolean prismatic_animation(fixed=false)
    "= true, if animation shall be enabled";
    parameter Real bodyCylinder1_body_I_2_1_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real relativePosition_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real prismatic_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_body_R_start_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Integer world_axisColor_z_3_(max=255, fixed=false, min=0)
    "Color of z-arrow";
    parameter Real bodyShape_body_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real der_bodyShape_body_w_a_1_(fixed=false, unit="rad/s2")
    "der(Absolute angular velocity of frame_a resolved in frame_a)";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles1_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles1_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real der_bodyCylinder1_body_w_a_1_(fixed=false, unit="rad/s2")
    "der(Absolute angular velocity of frame_a resolved in frame_a)";
    parameter Real relativePosition_zeroPosition_frame_resolve_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real revolute2_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_I_3_2_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real bodyCylinder_body_w_0_start_3_(fixed=false, unit="rad/s")
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real relativePosition_relativePosition_frame_resolve_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyShape_body_phi_2_(displayUnit="deg", unit="rad", fixed=false)
    "Dummy or 3 angles to rotate world frame into frame_a of body";
    parameter Real relativeVelocity_frame_a_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer relativeAngles1_sequence_2_(max=3, fixed=false, min=1)
    "Angles are returned to rotate frame_a around axes sequence[1], sequence[2] and finally sequence[3] into frame_b";
    parameter Modelica_SIunits_Distance rev_cylinderDiameter(fixed=false)
    "Diameter of cylinder representing the joint axis";
    parameter Integer bodyCylinder_sequence_start_2_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Integer bodyCylinder_body_cylinderColor_3_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real body_phi_dd_2_(fixed=false, unit="rad/s2") "= der(phi_d)";
    parameter Real relativeAngularVelocity1_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_I_1_2_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Modelica_SIunits_Length bodyCylinder1_length(fixed=false)
    "Length of cylinder";
    parameter Boolean bodyCylinder1_angles_fixed(fixed=false)
    "= true, if angles_start are used as initial values, else as guess values";
    parameter Real bodyCylinder_I_2_3_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real relativeAngularVelocity1_frame_b_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Boolean body_enforceStates(fixed=false)
    " = true, if absolute variables of body object shall be used as states (StateSelect.always)";
    parameter Real bodyCylinder1_frameTranslation_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Length world_gravityHeadWidth(fixed=false);
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer bodyCylinder1_body_sequence_start_1_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Integer bodyCylinder1_sequence_start_2_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngles_R_rel_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB relativeAngularVelocity1_relativeAngularVelocity_resolveInFrame(fixed=false)
    "Frame in which output vector w_rel is resolved (1: world, 2: frame_a, 3: frame_b, 4: frame_resolve)";
    parameter Real bodyCylinder1_body_Q_start_2_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Real bodyCylinder_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer body_sphereColor_3_(max=255, fixed=false, min=0)
    "Color of sphere";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Distance bodyCylinder_frameTranslation_width(fixed=false)
    " Width of shape";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frame_b_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_body_r_CM_1_(fixed=false, unit="m")
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Real world_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer rev_cylinderColor_3_(max=255, fixed=false, min=0)
    "Color of cylinder representing the joint axis";
    parameter Real relativePosition_zeroPosition_frame_resolve_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_R_rel_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Distance prismatic_boxWidth(fixed=false)
    "Width of prismatic joint box";
    parameter Real bodyShape_body_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Boolean bodyShape_body_useQuaternions(fixed=false)
    " = true, if quaternions shall be used as potential states otherwise use 3 angles as potential states";
    parameter Integer world_gravitySphereColor_2_(max=255, fixed=false, min=0)
    "Color of gravity sphere";
    parameter Real relativeVelocity_frame_b_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_r_0_3_(fixed=false, unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Boolean bodyCylinder_enforceStates(fixed=false)
    " = true, if absolute variables of body object shall be used as states (StateSelect.always)";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeAngularVelocity_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer world_axisColor_z_2_(max=255, fixed=false, min=0)
    "Color of z-arrow";
    parameter Real relativePosition_relativePosition_frame_resolve_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_Mechanics_MultiBody_Types_ShapeType bodyShape_frameTranslation_shapeType =  " Type of shape";
    parameter Real relativeVelocity_relativePosition_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean bodyCylinder1_w_0_fixed(fixed=false)
    "= true, if w_0_start are used as initial values, else as guess values";
    parameter Real bodyCylinder1_frameTranslation_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_zeroPosition_frame_resolve_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_body_phi_dd_1_(fixed=false, unit="rad/s2")
    "= der(phi_d)";
    parameter Real bodyShape_v_0_3_(fixed=false, unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real revolute2_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean revolute2_constantTorque_useSupport(fixed=false)
    "= true, if support flange enabled, otherwise implicitly grounded";
    parameter Real bodyShape_frameTranslation_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real der_bodyCylinder_body_frame_a_r_0_2_(fixed=false, unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    parameter Real der_body_w_a_2_(fixed=false, unit="rad/s2")
    "der(Absolute angular velocity of frame_a resolved in frame_a)";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean bodyCylinder_body_z_0_fixed(fixed=false)
    "= true, if z_0_start are used as initial values, else as guess values";
    parameter Real relativePosition_frame_b_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_phi_2_(displayUnit="deg", unit="rad", fixed=false)
    "Dummy or 3 angles to rotate world frame into frame_a of body";
    parameter Real der_bodyShape_body_v_0_3_(fixed=false, unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    parameter Real bodyCylinder_r_CM_2_(fixed=false, unit="m")
    "Position vector from frame_a to center of mass, resolved in frame_a";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativePosition_relativePosition_frame_resolve_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_I_1_2_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_R_rel_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real der_bodyCylinder_body_w_a_1_(fixed=false, unit="rad/s2")
    "der(Absolute angular velocity of frame_a resolved in frame_a)";
    parameter Real relativeVelocity_relativePosition_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB relativeVelocity_tansformRelativeVector_basicTransformVector_frame_r_out(fixed=false)
    "Frame in which vector r_out (= r_in in other frame) is resolved (1: world, 2: frame_a, 3: frame_b, 4: frame_resolve)";
    parameter Real relativeAngles1_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_frameTranslation_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_phi_1_(displayUnit="deg", unit="rad", fixed=false)
    "Dummy or 3 angles to rotate world frame into frame_a of body";
    parameter Integer relativeAngles1_sequence_3_(max=3, fixed=false, min=1)
    "Angles are returned to rotate frame_a around axes sequence[1], sequence[2] and finally sequence[3] into frame_b";
    parameter Real relativeVelocity_tansformRelativeVector_r_out_2_(fixed=false, unit="m/s")
    "Input vector r_in resolved in frame defined by frame_r_out";
    parameter Integer bodyCylinder_sequence_start_3_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Modelica_SIunits_Diameter bodyShape_body_sphereDiameter(fixed=false)
    "Diameter of sphere";
    parameter Real body_phi_dd_3_(fixed=false, unit="rad/s2") "= der(phi_d)";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_Q_4_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Real bodyCylinder_body_w_0_start_1_(fixed=false, unit="rad/s")
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB relativeVelocity_relativePosition_relativePosition_resolveInFrame(fixed=false)
    "Frame in which output vector r_rel is resolved (1: world, 2: frame_a, 3: frame_b, 4: frame_resolve)";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_w_rel_2_(fixed=false, unit="1/s")
    "Relative angular velocity vector between frame_a and frame_b resolved in frame defined by resolveInFrame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_R1_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_frameTranslation_widthDirection_3_(fixed=false, unit="1")
    " Vector in width direction of shape, resolved in frame_a";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_a_0_3_(fixed=false, unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_R_rel_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Modelica_SIunits_Distance bodyCylinder1_radius(fixed=false)
    "Radius of cylinder";
    parameter Real world_frame_b_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_resolve_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativePosition_relativePosition_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles_R_rel_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Integer bodyCylinder_sequence_angleStates_1_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real bodyCylinder1_body_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_r_CM_2_(fixed=false, unit="m")
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Real bodyCylinder_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real revolute2_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_w_a_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of frame_a resolved in frame_a";
    parameter Real bodyShape_body_Q_start_4_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Real prismatic_e_3_(fixed=false, unit="1")
    "Unit vector in direction of prismatic axis n";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Modelica_SIunits_Diameter bodyCylinder_body_cylinderDiameter(fixed=false)
    "Diameter of cylinder";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativePosition_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Length bodyShape_length(fixed=false)
    " Length of shape";
    parameter Real rev_R_rel_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder1_body_I_2_2_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real relativeVelocity_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real prismatic_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Length world_labelStart(fixed=false);
    parameter Real revolute2_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real prismatic_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real revolute2_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativePosition_relativePosition_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real der_bodyCylinder_body_frame_a_r_0_3_(fixed=false, unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    parameter Real der_body_w_a_1_(fixed=false, unit="rad/s2")
    "der(Absolute angular velocity of frame_a resolved in frame_a)";
    parameter Real bodyCylinder1_body_r_CM_3_(fixed=false, unit="m")
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Real bodyShape_frameTranslation_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real world_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_w_a_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of frame_a resolved in frame_a";
    parameter Real bodyCylinder1_body_phi_d_1_(fixed=false, unit="rad/s")
    "= der(phi)";
    parameter Real bodyCylinder1_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_b_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativePosition_relativePosition_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real world_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Integer revolute2_cylinderColor_3_(max=255, fixed=false, min=0)
    "Color of cylinder representing the joint axis";
    parameter Real relativePosition_relativePosition_frame_resolve_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_r_0_2_(fixed=false, unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    parameter Real bodyShape_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real der_bodyShape_body_v_0_2_(fixed=false, unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    parameter Real relativeAngles_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer world_ndim2(fixed=false);
    parameter Real rev_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_frameTranslation_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder1_body_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder1_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_b_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder1_body_z_0_start_2_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real relativePosition_relativePosition_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_R_rel_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean world_axisShowLabels(fixed=false)
    "= true, if labels shall be shown";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real world_frame_b_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_body_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_body_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_r_1_(fixed=false, unit="m")
    "Vector from frame_a to frame_b, resolved in frame_a";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_frameTranslation_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_I_3_1_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real der_bodyCylinder_body_w_a_2_(fixed=false, unit="rad/s2")
    "der(Absolute angular velocity of frame_a resolved in frame_a)";
    parameter Boolean bodyShape_body_z_0_fixed(fixed=false)
    "= true, if z_0_start are used as initial values, else as guess values";
    parameter Real bodyShape_body_w_0_start_1_(fixed=false, unit="rad/s")
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Position prismatic_s_offset(fixed=false)
    "Relative distance offset (distance between frame_a and frame_b = s_offset + s)";
    parameter Integer world_axisColor_x_1_(max=255, fixed=false, min=0)
    "Color of x-arrow";
    parameter Boolean world_animateGravity(fixed=false)
    "= true, if gravity field shall be visualized (acceleration vector or field center)";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Modelica_SIunits_Distance world_axisLength(fixed=false)
    "Length of world axes arrows";
    parameter Boolean bodyCylinder_useQuaternions(fixed=false)
    " = true, if quaternions shall be used as potential states otherwise use 3 angles as potential states";
    parameter Real bodyShape_body_R_start_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer world_axisColor_z_1_(max=255, fixed=false, min=0)
    "Color of z-arrow";
    parameter Real der_bodyShape_body_w_a_3_(fixed=false, unit="rad/s2")
    "der(Absolute angular velocity of frame_a resolved in frame_a)";
    parameter Real relativeVelocity_relativePosition_frame_b_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_frameTranslation_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_frameTranslation_r_shape_1_(fixed=false, unit="m")
    " Vector from frame_a to shape origin, resolved in frame_a";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Length world_defaultForceWidth(fixed=false)
    "Default for the fixed width of a shape represening a force (e.g., spring, bushing)";
    parameter Real bodyCylinder_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Integer bodyShape_sequence_angleStates_2_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real bodyCylinder_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real rev_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_SIunits_Mass bodyShape_body_m(fixed=false)
    "Mass of rigid body";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_R1_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean bodyShape_enforceStates(fixed=false)
    " = true, if absolute variables of body object shall be used as states (StateSelect.always)";
    parameter Real rev_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_resolve_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real prismatic_e_2_(fixed=false, unit="1")
    "Unit vector in direction of prismatic axis n";
    parameter Modelica_Mechanics_MultiBody_Types_ShapeType bodyCylinder1_frameTranslation_shapeType =  " Type of shape";
    parameter Real relativePosition_relativePosition_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_Mechanics_MultiBody_Types_ResolveInFrameAB relativeVelocity_tansformRelativeVector_basicTransformVector_frame_r_in(fixed=false)
    "Frame in which vector r_in is resolved (1: world, 2: frame_a, 3: frame_b, 4: frame_resolve)";
    parameter Real relativeVelocity_relativePosition_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_r_shape_3_(fixed=false, unit="m")
    " Vector from frame_a to shape origin, resolved in frame_a";
    parameter Real prismatic_frame_b_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_R_rel_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles1_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_frameTranslation_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean bodyShape_angles_fixed(fixed=false)
    "= true, if angles_start are used as initial values, else as guess values";
    parameter Real relativeVelocity_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Length length(fixed=false);
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_r_out_3_(fixed=false, unit="m/s")
    "Input vector r_in resolved in frame defined by frame_r_out";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder1_frameTranslation_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_body_R_start_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_r_CM_3_(fixed=false, unit="m")
    "Vector from frame_a to center of mass, resolved in frame_a";
    parameter Real bodyCylinder1_body_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real body_a_0_3_(fixed=false, unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    parameter Real prismatic_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_SIunits_Mass bodyCylinder_body_m(fixed=false)
    "Mass of rigid body";
    parameter Real relativeAngularVelocity_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean body_angles_fixed(fixed=false)
    "= true, if angles_start are used as initial values, else as guess values";
    parameter Real revolute2_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyCylinder_body_cylinderColor_2_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real prismatic_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Length world_defaultBodyDiameter(fixed=false)
    "Default for diameter of sphere representing the center of mass of a body";
    parameter Integer bodyCylinder1_sequence_angleStates_1_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Boolean bodyCylinder1_body_animation(fixed=false)
    "= true, if animation shall be enabled (show cylinder and sphere)";
    parameter Modelica_SIunits_Distance revolute2_cylinderDiameter(fixed=false)
    "Diameter of cylinder representing the joint axis";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Integer bodyCylinder1_sequence_angleStates_3_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real world_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_resolve_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyShape_body_phi_dd_1_(fixed=false, unit="rad/s2")
    "= der(phi_d)";
    parameter Real bodyCylinder1_body_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Modelica_Mechanics_MultiBody_Types_ShapeExtra bodyCylinder1_frameTranslation_extra(fixed=false, unit="1")
    " Additional parameter depending on shapeType (see docu of Visualizers.Advanced.Shape).";
    parameter Real relativePosition_relativePosition_frame_resolve_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativePosition_zeroPosition_frame_resolve_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Modelica_Mechanics_MultiBody_Types_GravityTypes world_gravityType(fixed=false)
    "Type of gravity field";
    parameter Real relativePosition_relativePosition_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real prismatic_frame_a_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real rev_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder1_body_R_start_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativePosition_relativePosition_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyShape_r_shape_1_(fixed=false, unit="m")
    " Vector from frame_a to shape origin, resolved in frame_a";
    parameter Real relativePosition_relativePosition_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder1_body_z_0_start_3_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_0 = der(w_0)";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Modelica_SIunits_Inertia bodyCylinder_body_I_11(fixed=false, min=0.0)
    " (1,1) element of inertia tensor";
    parameter Real relativeAngularVelocity1_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_r_2_(fixed=false, unit="m")
    "Vector from frame_a to frame_b, resolved in frame_a";
    parameter Real bodyCylinder_body_a_0_2_(fixed=false, unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    parameter Real relativePosition_zeroPosition_frame_resolve_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real revolute2_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Modelica_SIunits_Angle damper_phi_nominal(displayUnit="rad", fixed=false)
    "Nominal value of phi_rel (used for scaling)";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_w_rel_1_(fixed=false, unit="rad/s")
    "Relative angular velocity vector";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativePosition_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Diameter bodyCylinder1_body_sphereDiameter(fixed=false)
    "Diameter of sphere";
    parameter Real relativeAngularVelocity_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder1_body_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_R_start_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_frame_b_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real world_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder1_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real der_bodyShape_body_w_a_2_(fixed=false, unit="rad/s2")
    "der(Absolute angular velocity of frame_a resolved in frame_a)";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativePosition_relativePosition_frame_resolve_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativePosition_zeroPosition_frame_resolve_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real der_bodyCylinder1_body_v_0_3_(fixed=false, unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    parameter Integer bodyShape_sequence_angleStates_3_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real world_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_w_a_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of frame_a resolved in frame_a";
    parameter Real relativePosition_relativePosition_frame_resolve_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyCylinder1_body_sphereColor_2_(max=255, fixed=false, min=0)
    "Color of sphere";
    parameter Real relativePosition_frame_a_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean bodyCylinder1_body_z_0_fixed(fixed=false)
    "= true, if z_0_start are used as initial values, else as guess values";
    parameter Real relativePosition_relativePosition_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder1_frame_b_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real world_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_R1_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_widthDirection_1_(fixed=false, unit="1")
    " Vector in width direction of shape, resolved in frame_a";
    parameter Real bodyShape_r_0_3_(fixed=false, unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_frame_b_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_R_rel_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_r_shape_2_(fixed=false, unit="m")
    " Vector from frame_a to shape origin, resolved in frame_a";
    parameter Real prismatic_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_R_start_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyShape_frameTranslation_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_a_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder_body_z_a_start_2_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_a = der(w_a), i.e., time derivative of angular velocity resolved in frame_a";
    parameter Real relativeVelocity_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer bodyCylinder_sequence_start_1_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Modelica_SIunits_Diameter world_gravityArrowDiameter(fixed=false)
    "Diameter of gravity arrow";
    parameter Real bodyShape_a_0_3_(fixed=false, unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean damper1_useHeatPort(fixed=false)
    "=true, if heatPort is enabled";
    parameter Real bodyCylinder_body_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_R1_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_Q_start_1_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a at initial time";
    parameter Real bodyCylinder1_body_R_start_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyCylinder_body_sequence_angleStates_1_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_R_rel_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Integer bodyCylinder_sequence_angleStates_3_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_I_2_1_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real relativePosition_zeroPosition_frame_resolve_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real rev_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_body_phi_dd_2_(fixed=false, unit="rad/s2")
    "= der(phi_d)";
    parameter Real relativeAngles_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_resolve_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativePosition_relativePosition_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_frameTranslation_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real revolute2_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real body_phi_dd_1_(fixed=false, unit="rad/s2") "= der(phi_d)";
    parameter Integer bodyCylinder_color_3_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativePosition_relativePosition_frame_resolve_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngles_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer bodyShape_sequence_start_1_(max=3, fixed=false, min=1)
    "Sequence of rotations to rotate frame_a into frame_b at initial time";
    parameter Real bodyCylinder_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_b_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativePosition_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real der_body_frame_a_r_0_3_(fixed=false, unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    parameter Real bodyCylinder_frameTranslation_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativePosition_relativePosition_r_rel_2_(fixed=false, unit="m")
    "Relative position vector frame_b.r_0 - frame_a.r_0 resolved in frame defined by resolveInFrame";
    parameter Real bodyCylinder_body_g_0_3_(fixed=false, unit="m/s2")
    "Gravity acceleration resolved in world frame";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity_frame_a_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder_I_1_3_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_R_rel_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_frame_b_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyShape_body_R_start_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeVelocity_frame_b_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_Q_4_(fixed=false, unit="1")
    "Quaternion orientation object from world frame to frame_a (dummy value, if quaternions are not used as states)";
    parameter Real world_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_frame_a_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Modelica_SIunits_Inertia bodyCylinder_body_I_21(fixed=false, min=-1E+060)
    " (2,1) element of inertia tensor";
    parameter Modelica_SIunits_Inertia bodyCylinder_body_I_22(fixed=false, min=0.0)
    " (2,2) element of inertia tensor";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_SIunits_Length world_gravityHeadLength(fixed=false);
    parameter Real body_w_a_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of frame_a resolved in frame_a";
    parameter Real bodyCylinder_r_3_(fixed=false, unit="m")
    "Vector from frame_a to frame_b, resolved in frame_a";
    parameter Integer bodyCylinder1_body_sphereColor_3_(max=255, fixed=false, min=0)
    "Color of sphere";
    parameter Real bodyCylinder_body_z_a_start_1_(fixed=false, unit="rad/s2")
    "Initial values of angular acceleration z_a = der(w_a), i.e., time derivative of angular velocity resolved in frame_a";
    parameter Real prismatic_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_body_w_a_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of frame_a resolved in frame_a";
    parameter Real relativeAngularVelocity1_zeroPosition_frame_resolve_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real prismatic_frame_b_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_body_I_2_3_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real bodyShape_widthDirection_2_(fixed=false, unit="1")
    " Vector in width direction of shape, resolved in frame_a";
    parameter Boolean bodyShape_z_0_fixed(fixed=false)
    "= true, if z_0_start are used as initial values, else as guess values";
    parameter Real bodyCylinder1_body_angles_start_1_(displayUnit="deg", unit="rad", fixed=false)
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Real relativeAngles_R_rel_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_R_start_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder_body_I_1_1_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyShape_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngles1_R_rel_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativePosition_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativePosition_zeroPosition_frame_resolve_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngles1_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity1_frame_a_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real body_phi_1_(displayUnit="deg", unit="rad", fixed=false)
    "Dummy or 3 angles to rotate world frame into frame_a of body";
    parameter Real bodyCylinder1_frameTranslation_frame_b_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder_frameTranslation_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_a_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_tansformRelativeVector_frame_a_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Length world_headWidth(fixed=false);
    parameter Real bodyCylinder1_body_I_3_3_(fixed=false, unit="kg.m2")
    "inertia tensor";
    parameter Integer bodyCylinder_body_sphereColor_3_(max=255, fixed=false, min=0)
    "Color of sphere";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer world_axisColor_x_3_(max=255, fixed=false, min=0)
    "Color of x-arrow";
    parameter Boolean bodyCylinder_body_angles_fixed(fixed=false)
    "= true, if angles_start are used as initial values, else as guess values";
    parameter Real bodyShape_body_w_0_start_3_(fixed=false, unit="rad/s")
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Integer bodyCylinder_body_sequence_angleStates_2_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real relativeVelocity_relativePosition_zeroPosition_frame_resolve_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder_frameTranslation_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Modelica_SIunits_Mass bodyCylinder1_body_m(fixed=false)
    "Mass of rigid body";
    parameter Integer bodyCylinder_color_1_(max=255, fixed=false, min=0)
    "Color of cylinder";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_body_phi_start_2_(displayUnit="deg", unit="rad", fixed=false)
    "Potential angle states at initial time";
    parameter Real relativePosition_relativePosition_frame_resolve_R_T_2_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_body_R_start_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder1_body_R_start_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real revolute2_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativePosition_frame_b_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real bodyCylinder1_body_angles_start_3_(displayUnit="deg", unit="rad", fixed=false)
    "Initial values of angles to rotate frame_a around 'sequence_start' axes into frame_b";
    parameter Boolean bodyShape_body_angles_fixed(fixed=false)
    "= true, if angles_start are used as initial values, else as guess values";
    parameter Real relativeAngularVelocity_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real bodyCylinder1_r_CM_2_(fixed=false, unit="m")
    "Position vector from frame_a to center of mass, resolved in frame_a";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_R_T_3_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real prismatic_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_frame_a_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity1_frame_b_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeAngularVelocity1_frame_b_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real revolute2_R_rel_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativePosition_zeroPosition_frame_resolve_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_zeroPosition_frame_resolve_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_a_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_v_0_3_(fixed=false, unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real rev_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_zeroPosition_frame_resolve_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity1_relativeAngularVelocity_frame_resolve_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyCylinder1_body_w_0_start_3_(fixed=false, unit="rad/s")
    "Initial or guess values of angular velocity of frame_a resolved in world frame";
    parameter Real relativeAngles_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_r_rel_3_(fixed=false, unit="m")
    "Relative position vector frame_b.r_0 - frame_a.r_0 resolved in frame defined by resolveInFrame";
    parameter Real relativePosition_relativePosition_frame_b_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Integer world_axisColor_y_1_(max=255, fixed=false, min=0);
    parameter Modelica_SIunits_Length world_gravityLineLength(fixed=false);
    parameter Modelica_SIunits_Mass bodyCylinder_m(fixed=false)
    "Mass of cylinder";
    parameter Real bodyShape_body_frame_a_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real body_R_start_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_frame_a_t_2_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real body_phi_d_1_(fixed=false, unit="rad/s") "= der(phi)";
    parameter Real bodyCylinder1_I_1_2_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Boolean rev_animation(fixed=false)
    "= true, if animation shall be enabled (show axis as cylinder)";
    parameter Real bodyCylinder1_r_CM_1_(fixed=false, unit="m")
    "Position vector from frame_a to center of mass, resolved in frame_a";
    parameter Real bodyShape_a_0_2_(fixed=false, unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    parameter Real relativeAngles1_frame_b_r_0_3_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Modelica_SIunits_Diameter bodyShape_sphereDiameter(fixed=false)
    " Diameter of sphere";
    parameter Real bodyCylinder_I_1_2_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real bodyCylinder1_frame_b_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngularVelocity_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_R_rel_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Modelica_SIunits_Length world_gravityArrowLength(fixed=false)
    "Length of gravity arrow";
    parameter Real relativeVelocity_tansformRelativeVector_zeroPosition_frame_resolve_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real body_I_1_1_(fixed=false, unit="kg.m2") "inertia tensor";
    parameter Real bodyCylinder1_body_phi_3_(displayUnit="deg", unit="rad", fixed=false)
    "Dummy or 3 angles to rotate world frame into frame_a of body";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_R_T_1_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativePosition_relativePosition_frame_resolve_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_relativePosition_frame_a_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_frame_a_R_T_3_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real relativeVelocity_relativePosition_frame_a_f_2_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Modelica_Mechanics_MultiBody_Types_ShapeType bodyShape_shapeType =  " Type of shape";
    parameter Real relativePosition_relativePosition_frame_resolve_t_3_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_R_T_2_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyCylinder1_r_3_(fixed=false, unit="m")
    "Vector from frame_a to frame_b, resolved in frame_a";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real world_frame_b_R_T_1_1_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Boolean bodyCylinder1_enforceStates(fixed=false)
    " = true, if absolute variables of body object shall be used as states (StateSelect.always)";
    parameter Integer world_gravityArrowColor_2_(max=255, fixed=false, min=0)
    "Color of gravity arrow";
    parameter Real relativeAngularVelocity_relativeAngularVelocity_frame_resolve_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Integer bodyCylinder_sequence_angleStates_2_(max=3, fixed=false, min=1)
    " Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    parameter Real bodyCylinder_I_2_2_(fixed=false, unit="kg.m2")
    "Inertia tensor of cylinder with respect to center of mass, resolved in frame parallel to frame_a";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real bodyShape_body_phi_dd_3_(fixed=false, unit="rad/s2")
    "= der(phi_d)";
    parameter Real bodyCylinder_body_w_a_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of frame_a resolved in frame_a";
    parameter Real relativeVelocity_tansformRelativeVector_frame_b_f_1_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeAngularVelocity_frame_b_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real bodyCylinder_body_frame_a_R_T_2_3_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real prismatic_frame_a_R_w_2_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_resolve_t_1_(fixed=false, unit="N.m")
    "Cut-torque resolved in connector frame";
    parameter Real relativePosition_zeroPosition_frame_resolve_R_w_1_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativeAngles_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativePosition_relativePosition_frame_b_R_T_3_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real bodyShape_frame_a_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativePosition_relativePosition_frame_resolve_f_3_(fixed=false, unit="N")
    "Cut-force resolved in connector frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_b_R_T_1_2_(fixed=false)
    "Transformation matrix from world frame to local frame";
    parameter Real rev_frame_b_r_0_2_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    parameter Real relativeVelocity_relativePosition_relativePosition_frame_a_R_w_3_(fixed=false, unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    parameter Real relativePosition_zeroPosition_frame_resolve_r_0_1_(fixed=false, unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";

    Modelica_SIunits_Angle revolute2_phi
    "Relative rotation angle from frame_a to frame_b";
    Modelica_SIunits_AngularVelocity revolute2_w
    "First derivative of angle phi (relative angular velocity)";
    Modelica_SIunits_AngularVelocity rev_w
    "First derivative of angle phi (relative angular velocity)";
    Modelica_SIunits_Velocity prismatic_v(start=_prismatic_v_start)
    "First derivative of s (relative velocity)";
    Modelica_SIunits_Angle rev_phi(start=_rev_phi_start)
    "Relative rotation angle from frame_a to frame_b";
    Modelica_SIunits_Position prismatic_s(start=_prismatic_s_start)
    "Relative distance between frame_a and frame_b";

    Real der_revolute2_w(unit="rad/s2")
    "der(First derivative of angle phi (relative angular velocity))";
    Real der_prismatic_v(unit="m/s2")
    "der(First derivative of s (relative velocity))";
    Real der_prismatic_s(unit="m/s")
    "der(Relative distance between frame_a and frame_b)";
    Real der_rev_w(unit="rad/s2")
    "der(First derivative of angle phi (relative angular velocity))";
    Real der_revolute2_phi(unit="rad/s")
    "der(Relative rotation angle from frame_a to frame_b)";
    Real der_rev_phi(unit="rad/s")
    "der(Relative rotation angle from frame_a to frame_b)";

    Modelica_SIunits_Velocity damper1_v_rel "Relative velocity (= der(s_rel))";
    Real rev_frame_a_t_1_(unit="N.m") "Cut-torque resolved in connector frame";
    Real relativeAngles1_frame_b_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity1_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real rev_frame_b_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity1_frame_a_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyShape_frameTranslation_frame_b_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real relativePosition_relativePosition_r_rel_1_(unit="m")
    "Relative position vector frame_b.r_0 - frame_a.r_0 resolved in frame defined by resolveInFrame";
    Real der_bodyCylinder_body_v_0_1_(unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    Real bodyCylinder_frame_b_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real relativeAngularVelocity1_relativeAngularVelocity_w_rel_3_(unit="rad/s")
    "Relative angular velocity vector";
    Real der_body_frame_a_r_0_2_(unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    Real bodyCylinder1_frame_b_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real rev_R_rel_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real bodyShape_frame_a_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real bodyCylinder_frameTranslation_frame_a_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real relativeAngularVelocity_relativeAngularVelocity_frame_b_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frameTranslation_frame_b_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder1_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyShape_a_0_1_(unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    Real relativeAngles1_frame_b_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyCylinder1_frameTranslation_frame_b_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real body_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real prismatic_frame_a_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder_frameTranslation_frame_a_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_r_0_1_(unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    Real bodyCylinder_body_z_a_3_(unit="rad/s2")
    "Absolute angular acceleration of frame_a resolved in frame_a";
    Real bodyCylinder_body_frame_a_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real body_v_0_2_(unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    Real rev_frame_b_f_1_(unit="N") "Cut-force resolved in connector frame";
    Real relativeAngles_frame_b_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_frame_a_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_frameTranslation_frame_b_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real bodyCylinder_frame_b_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real relativePosition_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real body_frame_a_f_2_(unit="N") "Cut-force resolved in connector frame";
    Real relativeAngularVelocity1_w_rel_3_(unit="1/s")
    "Relative angular velocity vector between frame_a and frame_b resolved in frame defined by resolveInFrame";
    Modelica_SIunits_AngularVelocity damper_w_rel
    "Relative angular velocity (= der(phi_rel))";
    Real body_frame_a_f_1_(unit="N") "Cut-force resolved in connector frame";
    Real bodyCylinder_body_frame_a_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeVelocity_tansformRelativeVector_basicTransformVector_r_in_1_
    "Input vector resolved in frame defined by frame_r_in";
    Real bodyCylinder1_frameTranslation_frame_b_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real relativeAngularVelocity_relativeAngularVelocity_frame_b_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real der_bodyCylinder1_body_w_a_3_(unit="rad/s2")
    "der(Absolute angular velocity of frame_a resolved in frame_a)";
    Real bodyCylinder1_body_r_0_2_(unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    Modelica_SIunits_Position prismatic_support_s "Absolute position of flange";
    Real revolute2_R_rel_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_frameTranslation_frame_a_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity_relativeAngularVelocity_frame_b_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Modelica_SIunits_Angle damper_phi_rel(nominal=0.0001)
    "Relative rotation angle (= flange_b.phi - flange_a.phi)";
    Real relativeAngles1_R_rel_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real relativeAngles1_R_rel_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real revolute2_frame_b_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real bodyShape_frameTranslation_frame_b_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real revolute2_R_rel_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real relativeAngles1_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real revolute2_frame_b_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real rev_frame_b_f_2_(unit="N") "Cut-force resolved in connector frame";
    Real relativeAngularVelocity1_frame_b_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real bodyShape_frame_a_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeAngles_frame_b_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real revolute2_frame_b_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeVelocity_relativePosition_relativePosition_r_rel_1_(unit="m")
    "Relative position vector frame_b.r_0 - frame_a.r_0 resolved in frame defined by resolveInFrame";
    Real body_v_0_1_(unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    Modelica_SIunits_Torque revolute2_constantTorque_tau
    "Accelerating torque acting at flange (= -flange.tau)";
    Real bodyShape_frameTranslation_frame_a_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeAngles_R_rel_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real relativeAngularVelocity_relativeAngularVelocity_R_rel_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real bodyShape_v_0_1_(unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    Real relativeAngles1_R_rel_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real bodyCylinder_frame_b_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder_frameTranslation_frame_b_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real bodyCylinder_body_frame_a_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real bodyShape_frame_b_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Modelica_SIunits_AngularAcceleration damper_a_rel
    "Relative angular acceleration (= der(w_rel))";
    Real bodyCylinder_frame_a_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real relativeVelocity_relativePosition_r_rel_1_(unit="m")
    "Relative position vector resolved in frame defined by resolveInFrame";
    Real der_bodyCylinder_r_0_1_(unit="m/s")
    "der(Position vector from origin of world frame to origin of frame_a)";
    Real rev_frame_a_t_3_(unit="N.m") "Cut-torque resolved in connector frame";
    Real relativeAngles1_frame_b_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frameTranslation_frame_a_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity_relativeAngularVelocity_frame_b_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Modelica_Blocks_Interfaces_RealOutput const1_y
    "Connector of Real output signal";
    Real bodyShape_frameTranslation_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Modelica_SIunits_Torque rev_internalAxis_flange_tau
    "Cut torque in the flange";
    Real rev_frame_b_t_2_(unit="N.m") "Cut-torque resolved in connector frame";
    Real relativeAngles_R_rel_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real bodyShape_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real relativeAngularVelocity_frame_b_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real prismatic_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyShape_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyShape_frameTranslation_frame_a_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real bodyCylinder1_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyShape_frameTranslation_frame_a_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real revolute2_R_rel_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real bodyShape_frameTranslation_frame_b_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder_frame_b_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_a_0_2_(unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    Modelica_SIunits_Distance damper1_s_rel(nominal=0.0001)
    "Relative distance (= flange_b.s - flange_a.s)";
    Real bodyShape_body_v_0_1_(unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    Real prismatic_frame_a_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real rev_frame_b_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real bodyShape_frame_b_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real der_bodyCylinder_v_0_1_(unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    Real relativeVelocity_relativePosition_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyCylinder1_frameTranslation_frame_a_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real bodyShape_frame_a_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real der_bodyCylinder_body_w_a_3_(unit="rad/s2")
    "der(Absolute angular velocity of frame_a resolved in frame_a)";
    Real revolute2_frame_b_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frameTranslation_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyShape_body_r_0_1_(unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    Real bodyCylinder_body_frame_a_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real relativeAngularVelocity1_frame_a_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real revolute2_R_rel_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real der_relativeVelocity_relativePosition_relativePosition_frame_b_r_0_1_(unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    Real bodyCylinder_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyCylinder1_body_frame_a_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real body_z_a_3_(unit="rad/s2")
    "Absolute angular acceleration of frame_a resolved in frame_a";
    Real rev_frame_a_t_2_(unit="N.m") "Cut-torque resolved in connector frame";
    Real rev_frame_b_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real relativeAngles1_frame_b_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_Force prismatic_axis_f "Cut force directed into flange";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real bodyCylinder1_body_frame_a_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngles1_angles_1_(displayUnit="deg", unit="rad")
    "Angles to rotate frame_a into frame_b via 'sequence'";
    Real bodyCylinder1_body_frame_a_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real body_frame_a_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyShape_frameTranslation_frame_b_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder_frame_b_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngles1_R_rel_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real rev_frame_b_t_3_(unit="N.m") "Cut-torque resolved in connector frame";
    Real relativeAngularVelocity1_frame_b_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyCylinder_body_frame_a_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frameTranslation_frame_b_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frameTranslation_frame_a_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real bodyCylinder_frame_a_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real body_frame_a_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real bodyShape_frameTranslation_frame_a_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder1_frame_b_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeAngularVelocity_frame_b_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_a_0_1_(unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    Real relativeVelocity_tansformRelativeVector_basicTransformVector_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real prismatic_frame_a_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real revolute2_frame_a_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Modelica_SIunits_Torque damper_flange_b_tau "Cut torque in the flange";
    Real relativeAngles_angles_2_(displayUnit="deg", unit="rad")
    "Angles to rotate frame_a into frame_b via 'sequence'";
    Real bodyCylinder_frameTranslation_frame_b_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeAngularVelocity_frame_b_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real der_bodyCylinder1_v_0_1_(unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    Real relativeVelocity_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real relativeAngularVelocity_relativeAngularVelocity_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real revolute2_frame_a_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_frameTranslation_frame_b_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_Angle revolute2_constantTorque_phi
    "Angle of flange with respect to support (= flange.phi - support.phi)";
    Modelica_SIunits_Length force_s
    "Distance between flange and support (= flange.s - support.s)";
    Real relativeAngles_R_rel_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real bodyShape_body_frame_a_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Modelica_Blocks_Interfaces_RealInput force_f
    "Driving force as input signal";
    Real bodyCylinder_body_w_a_3_(unit="rad/s")
    "Absolute angular velocity of frame_a resolved in frame_a";
    Real bodyCylinder_frameTranslation_frame_b_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Modelica_SIunits_Torque revolute2_tau
    "Driving torque in direction of axis of rotation";
    Real bodyCylinder1_frame_b_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real bodyCylinder_frame_a_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyShape_frameTranslation_frame_b_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyShape_body_frame_a_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real revolute2_frame_a_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder_frameTranslation_frame_b_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder1_frame_a_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeAngles_angles_3_(displayUnit="deg", unit="rad")
    "Angles to rotate frame_a into frame_b via 'sequence'";
    Real bodyCylinder_body_frame_a_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frameTranslation_frame_b_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_body_a_0_1_(unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    Modelica_Mechanics_MultiBody_Types_SpecularCoefficient revolute2_specularCoefficient
    "Reflection of ambient light (= 0: light is completely absorbed)";
    Real relativeAngularVelocity1_frame_b_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real revolute2_frame_b_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real der_damper_w_rel(unit="rad/s2")
    "der(Relative angular velocity (= der(phi_rel)))";
    Real bodyCylinder_frameTranslation_frame_a_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_frameTranslation_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Modelica_SIunits_Angle rev_fixed_flange_phi
    "Absolute rotation angle of flange";
    Real bodyCylinder_frameTranslation_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Modelica_SIunits_Angle rev_angle "= phi_offset + phi";
    Real body_frame_a_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real der_bodyCylinder1_v_0_2_(unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    Real bodyShape_body_a_0_1_(unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    Modelica_Mechanics_MultiBody_Types_SpecularCoefficient bodyCylinder_body_specularCoefficient
    "Reflection of ambient light (= 0: light is completely absorbed)";
    Real bodyCylinder_frame_a_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real world_frame_b_f_2_(unit="N") "Cut-force resolved in connector frame";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real rev_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Modelica_Mechanics_MultiBody_Types_SpecularCoefficient rev_specularCoefficient
    "Reflection of ambient light (= 0: light is completely absorbed)";
    Real relativeAngles1_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real relativeAngles_R_rel_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real rev_frame_b_t_1_(unit="N.m") "Cut-torque resolved in connector frame";
    Real bodyShape_body_frame_a_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeVelocity_relativePosition_relativePosition_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyCylinder1_frame_b_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real der_relativeVelocity_der_r_rel_1__u
    "der(Connector of Real input signal)";
    Real bodyCylinder_frame_b_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_body_frame_a_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Modelica_SIunits_Force prismatic_fixed_flange_f
    "Cut force directed into flange";
    Modelica_SIunits_Position prismatic_fixed_flange_s
    "Absolute position of flange";
    Modelica_SIunits_Torque rev_axis_tau "Cut torque in the flange";
    Real bodyCylinder1_frameTranslation_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real revolute2_frame_a_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real relativeAngles1_angles_3_(displayUnit="deg", unit="rad")
    "Angles to rotate frame_a into frame_b via 'sequence'";
    Real bodyCylinder_body_frame_a_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real revolute2_frame_a_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frame_b_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Modelica_Blocks_Interfaces_RealOutput const_y
    "Connector of Real output signal";
    Real bodyCylinder_frameTranslation_frame_a_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real der_bodyCylinder1_r_0_1_(unit="m/s")
    "der(Position vector from origin of world frame to origin of frame_a)";
    Modelica_SIunits_Force damper1_flange_b_f "Cut force directed into flange";
    Real rev_frame_b_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_Position damper1_flange_b_s "Absolute position of flange";
    Real bodyCylinder1_v_0_2_(unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    Real bodyCylinder_frameTranslation_frame_b_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity1_relativeAngularVelocity_R_rel_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real relativeAngles1_frame_a_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frame_b_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real rev_R_rel_T_2_2_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_Torque rev_fixed_flange_tau "Cut torque in the flange";
    Real relativeAngularVelocity_relativeAngularVelocity_R_rel_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real relativeAngles1_frame_a_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_body_a_0_1_(unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    Real der_bodyShape_r_0_1_(unit="m/s")
    "der(Position vector from origin of world frame to origin of frame_a)";
    Modelica_SIunits_Force damper1_f "Forces between flanges (= flange_b.f)";
    Real relativeAngularVelocity1_frame_a_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_Angle damper_flange_a_phi
    "Absolute rotation angle of flange";
    Real body_frame_a_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_r_0_1_(unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    Modelica_SIunits_Angle revolute2_internalAxis_phi
    "External support angle (= flange.phi)";
    Real bodyCylinder1_frame_b_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder1_body_frame_a_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real relativeAngles1_frame_a_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_Angle damper_flange_b_phi
    "Absolute rotation angle of flange";
    Real bodyCylinder1_frameTranslation_frame_a_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativePosition_r_rel_1_(unit="m")
    "Relative position vector resolved in frame defined by resolveInFrame";
    Real bodyCylinder_frame_b_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real revolute2_frame_a_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real revolute2_frame_b_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_Force force_flange_f "Cut force directed into flange";
    Real relativeAngles_angles_1_(displayUnit="deg", unit="rad")
    "Angles to rotate frame_a into frame_b via 'sequence'";
    Real relativeVelocity_tansformRelativeVector_basicTransformVector_r_out_1_
    "Input vector r_in resolved in frame defined by frame_r_out";
    Real relativeAngularVelocity1_frame_b_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real rev_R_rel_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frameTranslation_frame_a_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Modelica_SIunits_Force prismatic_internalAxis_f
    "External support force (must be computed via force balance in model where InternalSupport is used; = flange.f)";
    Real body_frame_a_t_2_(unit="N.m") "Cut-torque resolved in connector frame";
    Modelica_SIunits_Position prismatic_internalAxis_s
    "External support position (= flange.s)";
    Real bodyCylinder1_frameTranslation_frame_b_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real relativeAngles1_angles_2_(displayUnit="deg", unit="rad")
    "Angles to rotate frame_a into frame_b via 'sequence'";
    Real relativeAngularVelocity_relativeAngularVelocity_frame_b_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real bodyCylinder1_frame_a_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_frame_b_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real der_damper1_v_rel(unit="m/s2") "der(Relative velocity (= der(s_rel)))";
    Modelica_SIunits_Angle revolute2_fixed_flange_phi
    "Absolute rotation angle of flange";
    Real relativeVelocity_v_rel_1_(unit="m/s")
    "Relative velocity vector resolved in frame defined by resolveInFrame";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngles1_frame_a_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity_relativeAngularVelocity_R_rel_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real relativeVelocity_tansformRelativeVector_r_in_1_
    "Input vector resolved in frame defined by frame_r_in";
    Real prismatic_frame_b_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Modelica_SIunits_Position prismatic_axis_s "Absolute position of flange";
    Real bodyCylinder1_body_v_0_2_(unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    Real body_frame_a_t_3_(unit="N.m") "Cut-torque resolved in connector frame";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity_relativeAngularVelocity_w_rel_3_(unit="rad/s")
    "Relative angular velocity vector";
    Real bodyCylinder1_frame_a_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeAngles_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyCylinder_frameTranslation_frame_a_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real bodyCylinder_frame_b_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeAngularVelocity1_frame_a_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_Angle rev_internalAxis_flange_phi
    "Absolute rotation angle of flange";
    Real bodyCylinder1_frameTranslation_frame_b_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real revolute2_frame_b_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyCylinder1_frameTranslation_frame_b_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder1_frame_b_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real rev_frame_a_f_2_(unit="N") "Cut-force resolved in connector frame";
    Real revolute2_R_rel_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real der_body_v_0_2_(unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    Real der_bodyCylinder1_body_frame_a_r_0_2_(unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    Real body_frame_a_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_Angle rev_axis_phi "Absolute rotation angle of flange";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_a_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_frame_a_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real der_bodyShape_body_v_0_1_(unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    Real bodyCylinder_body_r_0_1_(unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    Real world_frame_b_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeAngularVelocity1_relativeAngularVelocity_R_rel_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real bodyShape_frame_b_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeAngles1_frame_a_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real prismatic_frame_b_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real bodyCylinder1_body_z_a_3_(unit="rad/s2")
    "Absolute angular acceleration of frame_a resolved in frame_a";
    Real rev_frame_a_f_1_(unit="N") "Cut-force resolved in connector frame";
    Real bodyCylinder_frameTranslation_frame_b_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frameTranslation_frame_a_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity_relativeAngularVelocity_R_rel_T_1_1_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_Torque rev_tau
    "Driving torque in direction of axis of rotation";
    Real bodyCylinder1_frame_a_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real relativeAngularVelocity_w_rel_3_(unit="1/s")
    "Relative angular velocity vector between frame_a and frame_b resolved in frame defined by resolveInFrame";
    Real revolute2_frame_a_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real der_bodyCylinder1_body_frame_a_r_0_1_(unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    Modelica_Blocks_Interfaces_RealInput add_u1(unit="rad")
    "Connector of Real input signal 1";
    Modelica_Blocks_Interfaces_RealInput add_u2
    "Connector of Real input signal 2";
    Real rev_R_rel_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity1_frame_a_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real relativeAngles_frame_b_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real revolute2_frame_a_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Modelica_Blocks_Interfaces_RealOutput add_y
    "Connector of Real output signal";
    Real revolute2_frame_a_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real relativeAngularVelocity_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real revolute2_frame_b_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real der_body_w_a_3_(unit="rad/s2")
    "der(Absolute angular velocity of frame_a resolved in frame_a)";
    Real bodyCylinder1_frame_a_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real bodyCylinder1_frame_a_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_AngularAcceleration rev_a
    "Second derivative of angle phi (relative angular acceleration)";
    Real bodyCylinder1_body_a_0_2_(unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    Real bodyCylinder_body_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real der_body_v_0_1_(unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    Real bodyShape_frame_b_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Modelica_Mechanics_MultiBody_Types_SpecularCoefficient bodyShape_specularCoefficient
    "Reflection of ambient light (= 0: light is completely absorbed)";
    Real der_bodyShape_body_frame_a_r_0_1_(unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    Real bodyCylinder1_frame_b_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder1_body_w_a_3_(unit="rad/s")
    "Absolute angular velocity of frame_a resolved in frame_a";
    Real bodyCylinder1_frame_a_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Modelica_SIunits_Angle revolute2_constantTorque_flange_phi
    "Absolute rotation angle of flange";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_a_0_1_(unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    Modelica_SIunits_Torque revolute2_internalAxis_flange_tau
    "Cut torque in the flange";
    Real body_frame_a_t_1_(unit="N.m") "Cut-torque resolved in connector frame";
    Real relativeVelocity_tansformRelativeVector_r_out_1_(unit="m/s")
    "Input vector r_in resolved in frame defined by frame_r_out";
    Modelica_SIunits_Torque rev_support_tau "Cut torque in the flange";
    Real world_frame_b_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real body_a_0_1_(unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    Real bodyCylinder1_r_0_2_(unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    Real relativeAngularVelocity1_relativeAngularVelocity_R_rel_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real revolute2_frame_a_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real relativeAngularVelocity1_frame_a_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Modelica_Blocks_Interfaces_RealOutput add1_y
    "Connector of Real output signal";
    Real der_damper_phi_rel(unit="rad/s")
    "der(Relative rotation angle (= flange_b.phi - flange_a.phi))";
    Real bodyCylinder_frameTranslation_frame_b_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frameTranslation_frame_a_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frame_a_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_body_frame_a_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real rev_R_rel_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frame_b_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real der_damper1_s_rel(unit="m/s")
    "der(Relative distance (= flange_b.s - flange_a.s))";
    Modelica_Mechanics_MultiBody_Types_SpecularCoefficient prismatic_specularCoefficient
    "Reflection of ambient light (= 0: light is completely absorbed)";
    Real bodyCylinder_body_v_0_1_(unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    Modelica_SIunits_Angle rev_support_phi "Absolute rotation angle of flange";
    Real der_bodyCylinder1_body_v_0_1_(unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    Real relativeAngularVelocity1_relativeAngularVelocity_R_rel_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real bodyShape_r_0_1_(unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    Modelica_SIunits_Torque damper_flange_a_tau "Cut torque in the flange";
    Modelica_Blocks_Interfaces_RealInput add1_u1(unit="rad")
    "Connector of Real input signal 1";
    Modelica_Blocks_Interfaces_RealInput add1_u2
    "Connector of Real input signal 2";
    Real bodyCylinder1_frameTranslation_frame_b_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder_body_frame_a_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real relativeAngularVelocity_relativeAngularVelocity_R_rel_T_2_2_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_Angle revolute2_angle "= phi_offset + phi";
    Real body_frame_a_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Modelica_SIunits_Angle revolute2_internalAxis_flange_phi
    "Absolute rotation angle of flange";
    Real bodyCylinder_frameTranslation_frame_b_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real rev_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real relativeAngles_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real prismatic_frame_b_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real revolute2_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyCylinder_frame_b_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Modelica_SIunits_Position damper1_flange_a_s "Absolute position of flange";
    Real revolute2_frame_b_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real relativeAngularVelocity1_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real revolute2_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real body_r_0_2_(unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    Real der_body_frame_a_r_0_1_(unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    Real relativeAngularVelocity1_frame_b_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Modelica_Mechanics_MultiBody_Types_SpecularCoefficient bodyShape_body_specularCoefficient
    "Reflection of ambient light (= 0: light is completely absorbed)";
    Real bodyCylinder_frameTranslation_frame_a_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real body_a_0_2_(unit="m/s2")
    "Absolute acceleration of frame_a resolved in world frame (= der(v_0))";
    Modelica_Mechanics_MultiBody_Types_SpecularCoefficient bodyCylinder1_body_specularCoefficient
    "Reflection of ambient light (= 0: light is completely absorbed)";
    Real body_w_a_3_(unit="rad/s")
    "Absolute angular velocity of frame_a resolved in frame_a";
    Real world_frame_b_f_1_(unit="N") "Cut-force resolved in connector frame";
    Real relativeAngles_R_rel_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real prismatic_frame_a_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder_frameTranslation_frame_a_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyShape_frame_b_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Modelica_SIunits_AngularAcceleration revolute2_a
    "Second derivative of angle phi (relative angular acceleration)";
    Modelica_SIunits_Position force_flange_s "Absolute position of flange";
    Real bodyCylinder_v_0_1_(unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    Real bodyCylinder1_frame_b_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_Power damper_lossPower
    "Loss power leaving component via heatPort (> 0, if heat is flowing out of component)";
    Real bodyCylinder1_frameTranslation_frame_a_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_frame_b_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real relativeAngles_frame_b_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Modelica_Blocks_Interfaces_RealOutput relativeVelocity_der_r_rel_1__y
    "Connector of Real output signal";
    Modelica_SIunits_Force prismatic_support_f "Cut force directed into flange";
    Real world_frame_b_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real relativeAngularVelocity1_relativeAngularVelocity_R_rel_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_frame_a_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity_frame_b_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real relativeAngles_frame_b_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frame_a_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real prismatic_frame_b_t_2_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyShape_body_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real relativeVelocity_tansformRelativeVector_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real body_r_0_1_(unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    Modelica_SIunits_Power damper1_lossPower
    "Loss power leaving component via heatPort (> 0, if heat is flowing out of component)";
    Real relativeAngles1_R_rel_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real bodyShape_frameTranslation_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real revolute2_frame_b_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real relativeAngularVelocity1_frame_b_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_frameTranslation_frame_a_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder1_body_v_0_1_(unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    Real bodyCylinder1_frameTranslation_frame_b_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Real relativeAngles1_frame_b_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real prismatic_frame_a_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Modelica_SIunits_Torque revolute2_constantTorque_flange_tau
    "Cut torque in the flange";
    Real relativeAngles1_frame_a_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyCylinder1_body_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real prismatic_frame_b_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Modelica_SIunits_Angle rev_internalAxis_phi
    "External support angle (= flange.phi)";
    Real bodyCylinder1_frameTranslation_frame_a_R_T_2_2_
    "Transformation matrix from world frame to local frame";
    Real rev_frame_b_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real relativeAngularVelocity_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Modelica_SIunits_Force prismatic_f
    "Actuation force in direction of joint axis";
    Modelica_SIunits_Acceleration prismatic_a
    "Second derivative of s (relative acceleration)";
    Real der_bodyCylinder_body_frame_a_r_0_1_(unit="m/s")
    "der(Position vector from world frame to the connector frame origin, resolved in world frame)";
    Modelica_SIunits_Torque revolute2_internalAxis_tau
    "External support torque (must be computed via torque balance in model where InternalSupport is used; = flange.tau)";
    Real bodyCylinder_frame_a_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real bodyCylinder1_body_frame_a_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_body_r_0_1_(unit="m")
    "Position vector from origin of world frame to origin of frame_a";
    Modelica_SIunits_Torque damper_tau
    "Torque between flanges (= flange_b.tau)";
    Real relativePosition_relativePosition_frame_b_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyCylinder1_frameTranslation_frame_a_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real relativeAngularVelocity_frame_b_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_body_frame_a_f_1_(unit="N")
    "Cut-force resolved in connector frame";
    Modelica_SIunits_Force prismatic_internalAxis_flange_f
    "Cut force directed into flange";
    Real relativeAngles1_frame_b_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_frame_a_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real der_bodyCylinder1_body_v_0_2_(unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    Modelica_SIunits_Force damper1_flange_a_f "Cut force directed into flange";
    Real bodyShape_frame_a_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real revolute2_frame_a_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Modelica_Blocks_Interfaces_RealInput relativeVelocity_der_r_rel_1__u
    "Connector of Real input signal";
    Real bodyShape_frame_a_t_1_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder_frameTranslation_frame_a_t_3_(unit="N.m")
    "Cut-torque resolved in connector frame";
    Real bodyCylinder_frame_a_R_T_2_1_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder_frameTranslation_frame_a_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Modelica_Mechanics_MultiBody_Types_SpecularCoefficient bodyShape_frameTranslation_specularCoefficient
    "Reflection of ambient light (= 0: light is completely absorbed)";
    Real der_bodyShape_v_0_1_(unit="m/s2")
    "der(Absolute velocity of frame_a, resolved in world frame (= der(r_0)))";
    Real bodyCylinder1_frame_a_r_0_2_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Real bodyCylinder1_frameTranslation_frame_b_R_w_3_(unit="rad/s")
    "Absolute angular velocity of local frame, resolved in local frame";
    Real bodyShape_frameTranslation_frame_a_f_2_(unit="N")
    "Cut-force resolved in connector frame";
    Real relativeAngularVelocity_relativeAngularVelocity_frame_a_r_0_1_(unit="m")
    "Position vector from world frame to the connector frame origin, resolved in world frame";
    Modelica_SIunits_Torque rev_internalAxis_tau
    "External support torque (must be computed via torque balance in model where InternalSupport is used; = flange.tau)";
    Real bodyCylinder1_body_frame_a_R_T_1_1_
    "Transformation matrix from world frame to local frame";
    Modelica_SIunits_Position prismatic_internalAxis_flange_s
    "Absolute position of flange";
    Real der_bodyCylinder1_r_0_2_(unit="m/s")
    "der(Position vector from origin of world frame to origin of frame_a)";
    Real relativeAngularVelocity1_relativeAngularVelocity_frame_b_R_T_1_2_
    "Transformation matrix from world frame to local frame";
    Real bodyCylinder1_v_0_1_(unit="m/s")
    "Absolute velocity of frame_a, resolved in world frame (= der(r_0))";
    Real bodyCylinder_frameTranslation_frame_b_f_2_(unit="N")
    "Cut-force resolved in connector frame";

public
  parameter String fmi_instanceName="ModelInstance"
      annotation (Dialog(tab="FMI", group="Instance name"));
  parameter Boolean fmi_loggingOn=true
      annotation (Dialog(tab="FMI", group="Enable logging"));

protected
  constant Integer fmi_NumberOfContinuousStates=Functions.nx;
constant Integer fmi_NumberOfEventIndicators = 0;

  Functions.FmiInstance fmi=Functions.FmiInstance(fmi_instanceName, fmi_loggingOn);
    //Boolean fmi_z_positive[fmi_NumberOfEventIndicators];
  parameter Real fmi_AParamsAndStart(fixed=false);
    //Boolean fmi_StepEvent;

   // Real fmi_AStates;
    //parameter Real fmi_Initialized(fixed=false);
    //Integer fmi_DiscreteInputsSet = (if fmi_NewStates then 1 else 0);
    //Boolean fmi_DiscreteInputChanged=false;
    //Boolean fmi_NewStates;
    //Real fmi_TNext;

initial algorithm
  fmi_AParamsAndStart := 1;
 //Assign Parameters including Initial Values
// Start values

  //	1 Boolean parameters
  fmi_AParamsAndStart := Functions.fmiSetBoolean(
    fmi,
    {16777218},
    {world_driveTrainMechanics3D},
    1);

  //	282 Real parameters
  fmi_AParamsAndStart := Functions.fmiSetReal(
    fmi,
    {16777262,100665144,16777260,16777268,100665167,16777306,100665142,16777277,
      16777286,100663811,100664888,100664843,16777239,100663575,100663771,
      16777280,16777305,100665113,16777235,100665225,100664915,100664895,
      100663805,100664805,16777269,100663357,16777238,100664918,100663549,
      100664784,100665106,100663737,100663550,16777300,100663573,16777251,
      100665112,100663543,16777287,100664979,16777308,100665165,100665173,
      100663556,100663813,100663545,100665042,16777256,16777250,100663548,
      100663577,16777237,100664897,100663792,100665168,100664927,16777307,
      100665148,100663542,100665039,100665040,100663812,100663578,100663547,
      100665094,100663564,100663557,100663739,100663791,16777289,100664919,
      100664536,16777234,100663546,100663763,16777223,16777259,16777261,
      16777249,100665149,100663356,16777217,16777248,16777303,16777225,16777270,
      100663774,100664926,100663563,100663558,16777274,100664858,100665093,
      100663553,16777227,100664788,100665058,100663796,16777263,100665116,
      100665038,100663790,16777295,16777236,100664907,100663775,100663827,
      100664899,100663835,100664870,16777273,100664892,16777247,16777304,
      100665092,100663562,16777226,16777296,100664787,100663765,100664925,
      100663583,100663544,100665147,100663828,100665132,100663749,100663752,
      100663751,16777246,16777293,16777272,100663776,100665163,100664869,
      100663555,100665155,100663764,16777265,100663584,100663824,100663743,
      100665062,16777297,100663798,16777258,100664896,100663836,100663750,
      100663748,16777271,100663793,100664901,100665041,100664857,16777278,
      16777219,100663554,100664868,100663804,16777292,100664894,100663535,
      100664785,100665154,16777266,100665066,100665175,100663585,100664900,
      100663826,100663794,100663747,16777228,16777229,100663797,100664865,
      100663358,16777254,100664786,16777257,100663576,100665224,100663803,
      100663745,100663512,100664893,16777231,100665143,100664898,100663325,
      100663834,100663630,100664978,100663795,100664859,16777283,16777282,
      16777221,16777224,100663354,100663744,100663802,100664813,16777253,
      100665153,100664920,16777294,16777264,16777240,100664866,16777288,
      100663361,16777281,100665166,100664977,16777291,16777232,16777285,
      16777220,16777276,16777216,16777255,100664864,16777279,100665037,
      100663829,100665136,100665117,100663770,16777299,100663799,16777242,
      16777245,16777244,100663746,16777290,100663806,100664840,100665226,
      100665146,16777252,100665118,100664841,16777275,16777233,16777222,
      100663800,16777301,100663551,100664906,16777243,16777241,100663337,
      100664917,16777298,100663807,16777280,100664842,100663360,100664905,
      100665105,100665145,16777267,100664884,16777284,100663355,100665140,
      100663772,100665091,100665174,100663801,100665141,100665107,100665114,
      100663362,16777302,100664789,100663552,16777230,100663541,100664809},
    {bodyShape_I_31,bodyCylinder1_body_R_start_T_2_1_,bodyShape_I_33,
      bodyShape_w_0_start_2_,bodyCylinder1_body_Q_start_3_,
      bodyCylinder1_z_0_start_1_,bodyCylinder1_body_R_start_T_1_2_,add_k1,
      bodyCylinder_innerDiameter,bodyShape_body_z_a_start_1_,
      bodyCylinder_body_I_2_2_,bodyCylinder_body_I_33,body_r_CM_3_,
      body_Q_start_1_,bodyShape_body_w_0_start_2_,bodyCylinder1_frame_a_t_3_,
      bodyCylinder1_w_0_start_3_,bodyCylinder1_body_w_0_start_2_,rev_fixed_phi0,
      bodyCylinder1_frameTranslation_height,bodyCylinder_body_g_0_2_,
      bodyCylinder_body_R_start_T_1_3_,bodyShape_body_R_start_T_3_1_,
      bodyCylinder_I_1_1_,bodyShape_w_0_start_3_,world_lineWidth,body_r_CM_2_,
      bodyCylinder_body_Q_start_2_,body_I_3_3_,bodyCylinder_radius,
      bodyCylinder1_body_angles_start_2_,bodyShape_body_frame_a_f_2_,
      body_R_start_T_1_1_,bodyCylinder1_angles_start_1_,body_g_0_2_,
      body_w_0_start_3_,bodyCylinder1_body_w_0_start_1_,body_I_1_3_,
      bodyCylinder_density,bodyCylinder_frameTranslation_extra,
      bodyCylinder1_z_0_start_3_,bodyCylinder1_body_Q_start_1_,
      bodyCylinder1_body_phi_start_1_,body_R_start_T_3_1_,
      bodyShape_body_z_a_start_3_,body_I_2_2_,bodyCylinder1_m,bodyShape_r_CM_2_,
      body_w_0_start_2_,body_I_3_2_,body_Q_start_3_,body_r_CM_1_,
      bodyCylinder_body_R_start_T_2_2_,bodyShape_body_I_1_3_,
      bodyCylinder1_body_Q_start_4_,bodyCylinder_body_phi_start_3_,
      bodyCylinder1_z_0_start_2_,bodyCylinder1_body_R_start_T_3_2_,body_I_1_2_,
      bodyCylinder1_mo,bodyCylinder1_mi,bodyShape_body_z_a_start_2_,
      body_Q_start_4_,body_I_3_1_,bodyCylinder1_body_I_33,body_z_a_start_3_,
      body_R_start_T_3_2_,bodyShape_body_frame_a_t_1_,bodyShape_body_I_1_2_,
      bodyCylinder_angles_start_2_,bodyCylinder_body_Q_start_3_,
      revolute2_frame_b_t_3_,damper1_d,body_I_2_3_,
      bodyShape_body_angles_start_1_,world_defaultWidthFraction,bodyShape_I_22,
      bodyShape_I_21,body_w_0_start_1_,bodyCylinder1_body_R_start_T_3_3_,
      world_lineLength,world_mue,body_angles_start_3_,
      bodyCylinder1_w_0_start_1_,world_defaultSpecularCoefficient,
      bodyShape_z_0_start_1_,bodyShape_body_z_0_start_1_,
      bodyCylinder_body_phi_start_2_,body_z_a_start_2_,body_R_start_T_3_3_,
      bodyShape_frameTranslation_extra,bodyCylinder_body_angles_start_2_,
      bodyCylinder1_body_I_22,body_R_start_T_2_1_,world_defaultNm_to_m,
      bodyCylinder_I22,bodyCylinder1_I_1_1_,bodyShape_body_I_3_1_,
      bodyShape_I_32,bodyCylinder1_body_z_0_start_1_,bodyCylinder1_innerRadius,
      bodyShape_body_I_1_1_,bodyCylinder_z_0_start_2_,damper_d,
      bodyCylinder_body_z_a_start_3_,bodyShape_body_z_0_start_2_,
      bodyShape_body_Q_start_2_,bodyCylinder_body_R_start_T_3_1_,
      bodyShape_body_phi_start_2_,bodyCylinder_body_z_0_start_3_,
      bodyShape_extra,bodyCylinder_body_I_3_3_,body_angles_start_2_,
      bodyCylinder1_w_0_start_2_,bodyCylinder1_body_I_11,body_z_a_start_1_,
      world_defaultN_to_m,bodyCylinder_z_0_start_3_,bodyCylinder_mi,
      bodyShape_body_angles_start_3_,bodyCylinder_body_phi_start_1_,
      body_phi_start_1_,body_I_2_1_,bodyCylinder1_body_R_start_T_3_1_,
      bodyShape_body_Q_start_3_,bodyCylinder1_body_I_1_1_,bodyShape_body_I_33,
      bodyShape_body_I_32,bodyShape_body_I_31,body_angles_start_1_,
      bodyCylinder_w_0_start_3_,bodyShape_z_0_start_3_,
      bodyShape_body_z_0_start_3_,bodyCylinder1_body_g_0_2_,
      bodyCylinder_body_z_0_start_2_,body_R_start_T_2_3_,
      bodyCylinder1_body_z_a_start_3_,bodyShape_body_angles_start_2_,
      bodyShape_angles_start_2_,body_phi_start_2_,bodyShape_body_g_0_2_,
      bodyShape_body_r_CM_1_,bodyCylinder1_I_2_2_,bodyCylinder1_diameter,
      bodyShape_body_I_3_3_,bodyShape_I_11,bodyCylinder_body_R_start_T_2_1_,
      bodyShape_body_phi_start_3_,bodyShape_body_I_21,bodyShape_body_I_22,
      bodyShape_z_0_start_2_,bodyShape_body_I_2_1_,
      bodyCylinder_body_R_start_T_3_3_,bodyCylinder1_I22,
      bodyCylinder_body_angles_start_1_,add_k2,world_gravityArrowTail_1_,
      body_R_start_T_2_2_,bodyCylinder_body_z_0_start_1_,
      bodyShape_body_R_start_T_2_3_,bodyCylinder_w_0_start_2_,
      bodyCylinder_body_R_start_T_1_2_,body_specularCoefficient,
      bodyCylinder_innerRadius,bodyCylinder1_body_z_a_start_2_,
      bodyShape_angles_start_3_,bodyCylinder1_I_3_3_,
      bodyCylinder1_body_phi_start_3_,body_phi_start_3_,
      bodyCylinder_body_R_start_T_3_2_,bodyShape_body_Q_start_1_,
      bodyShape_body_I_2_2_,bodyShape_body_I_11,m_trolley,m_load,
      bodyShape_body_I_3_2_,bodyCylinder_body_w_0_start_2_,world_scaledLabel,
      body_z_0_start_3_,bodyCylinder_mo,bodyShape_r_CM_3_,body_Q_start_2_,
      bodyCylinder1_frameTranslation_width,bodyShape_body_R_start_T_2_2_,
      bodyShape_body_r_CM_3_,body_m,bodyCylinder_body_R_start_T_1_1_,w1_start,
      bodyCylinder1_body_R_start_T_1_3_,bodyCylinder_body_R_start_T_2_3_,
      world_axisDiameter,bodyShape_body_phi_start_1_,bodyShape_m,
      bodyCylinder_frameTranslation_height,bodyShape_body_I_2_3_,
      bodyCylinder_body_angles_start_3_,add1_k2,add1_k1,
      world_gravityArrowTail_3_,world_defaultFrameDiameterFraction,
      world_headLength,bodyShape_body_r_CM_2_,bodyShape_body_R_start_T_2_1_,
      bodyCylinder_I_3_3_,body_z_0_start_2_,bodyCylinder1_body_z_a_start_1_,
      bodyCylinder_body_Q_start_4_,bodyCylinder_z_0_start_1_,
      bodyShape_angles_start_1_,body_I_11,bodyCylinder_body_w_0_start_3_,
      bodyCylinder_angles_start_1_,world_gravityHeadWidth,
      relativeAngles1_guessAngle1,bodyCylinder1_body_Q_start_2_,
      bodyCylinder_frameTranslation_width,bodyCylinder_w_0_start_1_,w2_start,
      bodyCylinder_diameter,world_gravityArrowTail_2_,const_k,world_g,
      bodyShape_r_CM_1_,bodyCylinder_body_w_0_start_1_,revolute2_fixed_phi0,
      bodyCylinder1_radius,bodyShape_body_Q_start_4_,bodyCylinder1_body_I_2_2_,
      bodyCylinder1_body_z_0_start_2_,bodyShape_body_w_0_start_1_,
      bodyCylinder1_density,bodyShape_body_R_start_T_1_1_,body_I_33,body_I_32,
      body_I_31,bodyShape_body_m,bodyCylinder_angles_start_3_,
      bodyShape_body_R_start_T_3_2_,bodyCylinder_body_m,
      bodyCylinder1_frameTranslation_extra,bodyCylinder1_body_R_start_T_2_3_,
      body_z_0_start_1_,bodyCylinder1_body_z_0_start_3_,bodyCylinder_body_I_11,
      relativeAngles_guessAngle1,prismatic_fixed_s0,world_gravitySphereDiameter,
      bodyShape_body_R_start_T_1_2_,bodyCylinder1_angles_start_2_,
      body_R_start_T_1_2_,bodyCylinder_body_z_a_start_2_,body_I_21,body_I_22,
      world_gravityArrowDiameter,bodyCylinder_body_Q_start_1_,
      bodyCylinder1_innerDiameter,bodyShape_body_R_start_T_3_3_,
      revolute2_constantTorque_tau_constant,bodyCylinder_body_I_22,
      world_gravityHeadLength,bodyCylinder_body_z_a_start_1_,
      bodyCylinder1_body_angles_start_1_,bodyCylinder1_body_R_start_T_2_2_,
      bodyShape_w_0_start_1_,bodyCylinder_body_I_1_1_,const1_k,world_headWidth,
      bodyCylinder1_body_I_3_3_,bodyShape_body_w_0_start_3_,
      bodyCylinder1_body_m,bodyCylinder1_body_phi_start_2_,
      bodyShape_body_R_start_T_1_3_,bodyCylinder1_body_R_start_T_1_1_,
      bodyCylinder1_body_angles_start_3_,bodyCylinder1_body_w_0_start_3_,
      world_gravityLineLength,bodyCylinder1_angles_start_3_,bodyCylinder_m,
      body_R_start_T_1_3_,phi2_start,body_I_1_1_,bodyCylinder_I_2_2_},
    1);
  //	4 Real start parameters

  fmi_AParamsAndStart := Functions.fmiSetReal(
    fmi,
    {352321536,33554434,33554433,33554432},
    {_u_start,_rev_phi_start,_prismatic_v_start,_prismatic_s_start},
    1);

  Functions.fmiInitialize(fmi, 1);

initial equation

equation

  annotation (Icon(
      Rectangle(extent=[-100,100; 100,-100], style(
          color=10,
          rgbcolor= { 95,95,95},
          thickness=2)),
      Text(
        extent=[-100,20; 100,-20],
        string="%name",
        style(
          color=3,
          rgbcolor= { 0,0,255},
          thickness=2,
          fillPattern=1)),
      Text(
        extent=[-80,-72; 80,-94],
        style(
          color=10,
          rgbcolor= { 95,95,95},
          thickness=2,
          fillPattern=1),
        string="FMI import")),
        uses(Modelica(version="3.2")));
end DoublePendulum;
