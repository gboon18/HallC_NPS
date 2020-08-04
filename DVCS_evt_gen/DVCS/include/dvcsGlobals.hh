
struct dvcsGlobals
{
  static int    run_number;
  static double Ebeam; // Beam energy in GeV's
  static double HMS_angle; // central angle of HMS window
  static double HMS_momentum; // HMS momentum central value
  static double Calo_distance; // Calorimeter distance
  static double Calo_angle;    // Calorimeter angle
  static int    target_type;   // target type p or deuteron
  static int    target_gen_proc_type; // target type for determining process (like DVCS on p = 0, n = 1, or d = 2)
  static double target_density; // density of target
  static double target_offset; // shift of the target center from the SC center
  static double target_length; // target lengt
  static bool   hit_HMS_CALO_flag;

  static bool   SM_field_flag; //Sweeping Magnet field ON||OFF
  static double SM_field_str; //Sweeping Magnet field strength

  static double SM_angle; //Sweeper Magnet angle
};
