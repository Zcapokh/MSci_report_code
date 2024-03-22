/*! @file Force_drag_tiegcm.cpp
        @author Santosh Bhattarai
        @date 30 January 2020
        @brief UCL ODL implementation of drag force model using TIE-GCM
 */

#include "../include/Force_drag_tiegcm.h"


// This is the file in ODL that is an adapted version of the one written by Alex Forsyth
// Key changes are the the parameters, dates, time units, and filepaths


void Force_drag_tiegcm::setup(const Resident_constants &rso_const,
                              std::shared_ptr<Resident_variables> in_state) {
  std::cout << "Force_drag_tiegcm::setup [1] \n";
  Force_drag::setup(rso_const, in_state);
  loadfiles();
}

void Force_drag_tiegcm::compute_acceleration() {
  if (state->geodetic.alt < 250.0 || state->geodetic.alt > 550.0) {
    std::stringstream error;
    error.precision(16);
    error << "Force_drag: Altitude too low or high, " << state->geodetic.alt
          << " km.";
    state->errors.push_back(error.str());
  }

  if (state->eci.epoch.get_MJD_UTC() < 54640.0 ||
      state->eci.epoch.get_MJD_UTC() > 54643) {
    std::stringstream error;
    error.precision(16);
    error << "Force_drag: Time not within the bounds of the data, "
          << state->eci.epoch.get_MJD_UTC() << " : Time";
    state->errors.push_back(error.str());
  }

  if ((state->geodetic.lat) * 180 / M_PI < -87.5 ||//-58.76 ||
      (state->geodetic.lat) * 180 / M_PI > 87.5){//56.24) {
    std::stringstream error;
    error.precision(16);
    error << "Force_drag: Latitude not within the bounds of the data, "
          << state->geodetic.lat << " degrees";
    state->errors.push_back(error.str());
  }

  std::cout << "Force_drag_tiegcm::compute_acceleration [2] \n";
  double rho = get_density(state->eci.epoch.get_MJD_UTC(),
                           (state->geodetic.alt) * 100000,
                           (state->geodetic.lat) * 180 / M_PI,
                           (state->geodetic.lon) * 180 / M_PI) *
               1E3;
  state->atmos_density = rho;
  std::cout << "Time: \n";
  state->eci.epoch.print_MJD_UTC();
  std::cout << "Altitude (km): " << state->geodetic.alt << std::endl;
  std::cout << "Latitude: " << state->geodetic.lat << std::endl;
  std::cout << "Longitude: " << state->geodetic.lon << std::endl;
  std::cout << "Density: " << rho << "kg/m^3" << std::endl;

  a_ecef = minus500C_dAm * rho * state->ecef_v * state->ecef_rso_vel;
  state->total_a_ecef += a_ecef;
}

void Force_drag_tiegcm::loadfiles() {
  std::string Densityfilepath = "../res/dietrich/denarray_24_control.txt";
  std::ifstream in_densityfile(Densityfilepath);
  if (in_densityfile.is_open()) {
    double X;
    std::string line;
    std::string row;
    std::string element;
    while (getline(in_densityfile, line)) {
      std::stringstream row(line);
      //std::cout<<line<<std::endl;
      std::vector<double> vectorelement;
      while (getline(row, element, ',')) {
        std::istringstream o(element);
        o >> X;
        vectorelement.push_back(X);
     }
      densityfile.push_back(vectorelement);
    }
    in_densityfile.close();
  }

  std::string Heightfilepath = "../res/dietrich/altarray_24_control.txt";
  std::ifstream in_heightfile(Heightfilepath);
  if (in_heightfile.is_open()) {
    double X;
    std::string line;
    std::string row;
    std::string element;
    while (getline(in_heightfile, line)) {
      std::stringstream row(line);
      std::vector<double> vectorelement;
      while (getline(row, element, ',')) {
        std::istringstream o(element);
        o >> X;
        vectorelement.push_back(X);
      }
      heightfile.push_back(vectorelement);
    }
    in_heightfile.close();
  }
  std::cout << "Force_drag_tiegcm::loadfiles [3] \n";
}

double Force_drag_tiegcm::get_density(double time, double altitude,
                                      double latitude, double longitude) {
  std::vector<double> V;
  //// identifies the lowest index limits of longitude, latitude and time for
  /// the variable - the first element is given by 1
  int lo1 = ((floor(longitude * 0.2) / 0.2) + 180) / 5 ;
  int la1 = ((floor((latitude + 3.75) * 0.2) / 0.2) + 55) / 5 + 1;
  //int la1 = static_cast<int>(((floor((latitude + 3.75) * 0.2) / 0.2) + 55) / 5 + 1);

  double time1 = floor(((time - 54640) * 24 * 60) / 60) + 1;
  //double time1 = floor(((time - 54640.08839299) * 24 * 60) / 60) + 1;
  std::cout << "lo1 is " << lo1 << ", la1 is " << la1 << " and time1 is: " <<
  time1 << " " << std::endl; std::cout<< longitude << " :lo1: " << lo1 << ", " <<
  latitude <<" :la1: " <<la1 <<", "<< altitude << " :altitude: " << altitude
  <<", "<< time << " :time: " << time1 << std::endl; std::cout <<
  densityfile[lo1][la1]
  << std::endl;
  if (time1 < 0) {
    time1 = 1;
  }

  for (double i = time1; i <= time1 + 1; i++) {
    // runs through each of the two closest longitude points
    for (double lo = lo1; lo <= lo1 + 1; lo++) {
      if (lo == 72) {
        lo = 1;
      }
      // runs through each of the two closest latitude points
      for (double la = la1; la <= la1 + 1; la++) {
        double ln = 216 * (i - 1) + (la - 1) * 9;

        // int lo1_int =  static_cast<int>(ln);
        // int ln_int = static_cast<int>(ln);
        // cout << ln << endl;
        // cout << lo << endl;
        // cout << Density << endl;
        // cout << Height << endl;
        double Height = heightfile[ln][lo];
        double Density = densityfile[ln][lo];
        // cout << Height << endl;
        double lowheight;
        double highheight;
        double lowdensity;
        double highdensity;

        // sees if the extracted value is less than the height we are looking
        // at. If so store it as the current lowest height value and the
        // corrisponding density value
        if (Height < altitude) {
          lowheight = Height;
          lowdensity = Density;
          // cout << "PASS" << endl;
        }
        // if the value extracted is greather than the height we are looking at
        // the height value is not contained within the height database
        else if (Height >= altitude) {
          std::cout << "ERROR: the height value: " << altitude
                    << "km is not contained within the database "
                       "'SkimedHeights.txt', check raw data around row number: "
                    << ln << " and column number: " << lo1 << std::endl;
        }

        for (int j = 1; j < 8; j++) {
          Height = heightfile[ln + j][lo];
          // cout << Height << endl;
          Density = densityfile[ln + j][lo];
          if (Height < altitude) {
            lowheight = Height;
            lowdensity = Density;
          } else if (Height >= altitude) {
            highheight = Height;
            highdensity = Density;

            // cout << lowheight << " & " << highheight << endl;
            // cout << lowdensity << " & " << highdensity << endl;
            double vectorelement = exp(
                (altitude - lowheight) * (log(highdensity) - log(lowdensity)) /
                    (highheight - lowheight) +
                log(lowdensity));
            V.push_back(vectorelement);
            // cout << vectorelement << endl;
            break;
          }
        }
      }
      if (lo == 1 && longitude > 175.0) {
        lo = 72;
      }
    }
  }
  double DENSITY = interpolate(latitude, longitude, time, V);
  return DENSITY;
}

double Force_drag_tiegcm::interpolate(double latitude, double longitude,
                                      double time, std::vector<double> v) {
  double longitudelower = floor(longitude * 0.2) / 0.2;
  double latitudelower = floor((latitude + 3.75) * 0.2) / 0.2 - 3.75;
  // cout << longitudelower << "   " << latitudelower << endl;
  double longxlower = (v[1] - v[0]) * (latitude - latitudelower) / 5 + v[0];
  double latxlower = (v[3] - v[1]) * (longitude - longitudelower) / 5 + v[1];
  double latxupper = (v[2] - v[0]) * (longitude - longitudelower) / 5 + v[0];
  double longxupper = (v[3] - v[2]) * (latitude - latitudelower) / 5 + v[2];
  double realtime = ((time - 54640) * 24 * 60) / 60;
  double time1 = floor(((time - 54640) * 24 * 60) / 60);
  // cout << longitudelower << "   " << latitudelower << endl;
  double longxlower2 = (v[5] - v[4]) * (latitude - latitudelower) / 5 + v[4];
  double latxlower2 = (v[7] - v[5]) * (longitude - longitudelower) / 5 + v[5];
  double latxupper2 = (v[6] - v[4]) * (longitude - longitudelower) / 5 + v[4];
  double longxupper2 = (v[7] - v[6]) * (latitude - latitudelower) / 5 + v[6];

  double massdensity1 =
      (((latxlower - latxupper) * (latitude - latitudelower) / 5 + latxupper) +
       ((longxupper - longxlower) * (longitude - longitudelower) / 5 +
        longxlower)) /
      2;

  double massdensity2 =
      (((latxlower2 - latxupper2) * (latitude - latitudelower) / 5 +
        latxupper2) +
       ((longxupper2 - longxlower2) * (longitude - longitudelower) / 5 +
        longxlower2)) /
      2;

  double massdensityx =
      (massdensity2 - massdensity1) * (realtime - time1) + massdensity1;

  return massdensityx;
  std::cout << "Force_drag_tiegcm::interpolation [5] \n";
}
