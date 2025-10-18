/**
 *  ------------------------------------------------
Calculation of communication latency between a
ground station and a satellite.
------------------------------------------------
 */

#include <iostream>
#include <string>
#include <cctype>
#include <cmath>
#include <chrono>
#include <vector>

const int SPEED_OF_LIGHT = 299792458; //meters per second
const int EARTH_RADIUS = 6371000; //meters
const double MU = 3.986004418e14; //Earth's gravitational parameter (m^3 / s^2)
const double EARTH_ROTATIONS_PER_DAY = 1.0027379093;   // Number of rotations the earth makes in a day relative to the stars
const double e = 2.71828182845904523536;
const double PI = 3.14159265358979323846;

struct SatelliteMetadata {
    int satellite_catalog_number;
    char classification;
    int international_designator_year;
    int international_designator_launch_number;
    std::string international_designator_piece;
    int epoch_year;
    double epoch_day;
    double drag_rate;
    double drag_acceleration;
    double drag_term;
    int ephemeris_type;
    int element_number;
};

struct OrbitalParams {
    int satellite_catalog_number;
    double inclination;
    double right_ascension;
    double eccentricity;
    double argument_of_perigee;
    double mean_anomaly;
    double mean_motion;
    int revolutions;
};

//Function headers
SatelliteMetadata parseTLEline1(const std::string &TLE_line1);
OrbitalParams parseTLEline2(const std::string &TLE_line2);
double convertRadians(const double degrees);
double convertDegrees(const double radians);
double parseTLESci(const std::string &num);
void display_data(const std::string &name, const SatelliteMetadata &metadata, const OrbitalParams &params, const double longitude, 
    const double latitude);
double calculate_true_anomaly(const double mean_motion, const double mean_anomaly, const int epoch_year, const double epoch_day, const double eccentricity);
double calculate_eccentric_anomaly(const double mean_motion_final, const double eccentricity);
std::vector<std::vector<double>> calculate_satellite_position_vector(const double inclination, const double right_ascension, 
    const double argument_of_perigee, const double mean_motion, const double true_anomaly, const double eccentricity);
std::vector<std::vector<double>> multiply_two_matrices(const std::vector<std::vector<double>> &matrix1, const std::vector<std::vector<double>> &matrix2);
std::vector<std::vector<double>> calculate_ground_position_vector(const double latitude, const double longitude);
double calculate_latency_distance(const std::vector<std::vector<double>> &satellite_position_vector, const std::vector<std::vector<double>> &ground_position_vector);
void display_results(const double true_anomaly, const double latency_distance, const double latency);

SatelliteMetadata parseTLEline1(const std::string &TLE_line1){
    /**
     * Parses data from first line of TLE
     */
    SatelliteMetadata data;
    data.satellite_catalog_number = std::stoi(TLE_line1.substr(2, 5));
    data.classification = TLE_line1[7];
    data.international_designator_year = std::stoi(TLE_line1.substr(9,2));
    data.international_designator_launch_number = std::stoi(TLE_line1.substr(11,3));
    data.international_designator_piece = TLE_line1.substr(14,3);
    data.epoch_year = std::stoi(TLE_line1.substr(18,2));
    data.epoch_day = std::stod(TLE_line1.substr(20, 12));
    data.drag_rate = parseTLESci(TLE_line1.substr(33,10));
    data.drag_acceleration = parseTLESci(TLE_line1.substr(44, 8));
    data.drag_term = parseTLESci(TLE_line1.substr(53, 8));
    data.ephemeris_type = TLE_line1[62];
    data.element_number = std::stoi(TLE_line1.substr(64, 4));
    return data;
}

OrbitalParams parseTLEline2(const std::string &TLE_line2){
    /**
     * Parses data from second line of TLE
     */
    OrbitalParams params;
    params.satellite_catalog_number = std::stoi(TLE_line2.substr(2, 5));
    params.inclination = std::stod(TLE_line2.substr(8, 8));
    params.right_ascension = std::stod(TLE_line2.substr(17, 8));
    params.eccentricity = std::stod("0." + TLE_line2.substr(26, 7));
    params.argument_of_perigee = std::stod(TLE_line2.substr(34, 8));
    params.mean_anomaly = std::stod(TLE_line2.substr(43, 8));
    params.mean_motion = std::stod(TLE_line2.substr(52,11));
    params.revolutions = std::stoi(TLE_line2.substr(63, 5));

    return params;

}

double convertRadians(const double degrees){
    /*Converts given angle in degrees to radians*/
    return degrees * (PI/180);
}

double convertDegrees(const double radians){
    /*Converts given angle in radians to degrees*/
    return radians * (180/PI);
}

double parseTLESci(const std::string &num){
    /**
     * Handles specific formatting instructions for certain TLE
     * data such as assuming a decimal point in the beginning and
     * +/- signs indicating scientific notation.
     */
    std::string parsed_num = num;
    size_t insert_index = (num[0] == '+' || num[0] == '-') ? 1 : 0;
    if (num[1] != '.'){
        parsed_num.insert(insert_index, "0.");
    }
    size_t last = parsed_num.find_last_of("+-");
    if (last != std::string::npos && last > 0 && last != 1){
        parsed_num.insert(last, "E");
    }

    return std::stod(parsed_num);

}

void display_data(const std::string &name, const SatelliteMetadata &metadata, const OrbitalParams &params, const double longitude, const double latitude){
    /**
     * Displays inputted satelite and ground station data.
     */
    std::cout << std::endl << "---------------------------------------" << std::endl;
    std::cout << "Satellite name: " << name << std::endl << std::endl;
    std::cout << "Inclination: " << params.inclination << std::endl;
    std::cout << "Right Ascension: " << params.right_ascension << std::endl;
    std::cout << "Eccentricity: " << params.eccentricity << std::endl;
    std::cout << "Argument of Perigee: " << params.argument_of_perigee << std::endl;
    std::cout << "Mean Anomaly: " << params.mean_anomaly << std::endl;
    std::cout << "Mean Motion: " << params.mean_motion << std::endl;
    std::cout << "Revolutions: " << params.revolutions << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "Latitude: " << latitude << std::endl;
    std::cout << "Longitude: " << longitude << std::endl;
    std::cout << "---------------------------------------" << std::endl;
}

double calculate_true_anomaly(const double mean_motion, const double mean_anomaly, const int epoch_year, const double epoch_day, const double eccentricity){
    /**
     * Calculates the satellites true anomaly angle.
     */
    double mean_motion_seconds = mean_motion / 86400; //  convert to revs/s
    int full_year = (epoch_year > 56)  ?  (epoch_year + 1900)  : (epoch_year + 2000);
    double epoch_satellite_time = ((full_year - 1970) * 31536000) + (epoch_day * 86400);
    std::chrono::system_clock::time_point tp = std::chrono::system_clock::now();
    std::chrono::system_clock::duration dtn = tp.time_since_epoch();
    double current_epoch = dtn.count() * std::chrono::system_clock::period::num / std::chrono::system_clock::period::den;

    double delta_t = current_epoch - epoch_satellite_time;

    double mean_motion_final = std::fmod( mean_anomaly + (mean_motion_seconds * delta_t), 360);

    double eccentric_anomaly = calculate_eccentric_anomaly(mean_motion_final, eccentricity);

    double true_anomaly_rad = std::acos((eccentricity - std::cos(convertRadians(eccentric_anomaly)))/(eccentricity * std::cos(convertRadians(eccentric_anomaly)) - 1));
    double true_anomaly_deg = convertDegrees(true_anomaly_rad);

    return true_anomaly_deg;
}

double calculate_eccentric_anomaly(const double mean_motion_final, const double eccentricity){
    /**
     * Calculates the eccentric anomaly by Solving equation M = E - esin(E) for E using 
     * the Newton-Rhapson method.
     */

    double E = mean_motion_final; // initial guess
    double tol = 1e-8;
    int count {};
    double dE {};
    do {
        count++;
        double E_rad = convertRadians(E);
        double dE = (E - eccentricity*std::sin(E_rad) - mean_motion_final) / (1 - eccentricity*std::cos(E_rad));
        E = E - dE;
    } while (count < 10000 && std::abs(dE) > tol);

    return E;
}

std::vector<std::vector<double>> calculate_satellite_position_vector(const double inclination, const double right_ascension, 
    const double argument_of_perigee, const double mean_motion, const double true_anomaly, const double eccentricity){
    /**
     * Calculates position vector in meters from center of earth to satellite's current position.
     */
    double inclination_rad = convertRadians(inclination);
    double right_ascension_rad = convertRadians(right_ascension);
    double argument_of_perigee_rad = convertRadians(argument_of_perigee);

    std::vector<std::vector<double>> Rkright_ascension_matrix = {
        {std::cos(right_ascension_rad), -std::cos(right_ascension_rad), 0},
        {std::sin(right_ascension_rad), std::cos(right_ascension_rad), 0},
        {0,0,1}
    };

    std::vector<std::vector<double>> RIinclination_matrix = {
        {1,0,0},
        {0, std::cos(inclination_rad), -std::sin(inclination_rad)}, 
        {0 , std::sin(inclination_rad), std::cos(inclination_rad)}
    };

    std::vector<std::vector<double>> Rkargument_of_perigee_matrix = {
        {std::cos(argument_of_perigee_rad), -std::sin(argument_of_perigee_rad), 0},
        {std::sin(argument_of_perigee_rad), std::cos(argument_of_perigee_rad), 0},
        {0,0,1}
    };

    std::vector<std::vector<double>> intermediate_result = multiply_two_matrices(Rkright_ascension_matrix, RIinclination_matrix);
    std::vector<std::vector<double>> transformation_matrix = multiply_two_matrices(intermediate_result, Rkargument_of_perigee_matrix);

    double mean_motion_rad_seconds = (mean_motion * 2 * PI)/ 86400;
    double semi_major_axis = std::pow(MU / (mean_motion_rad_seconds * mean_motion_rad_seconds), 1.0/3.0);

    double radius_mag = (semi_major_axis * (1 - (eccentricity * eccentricity))) / (1 + (eccentricity * std::cos(convertRadians(true_anomaly))));

    std::vector<std::vector<double>> radius_ijk = {
        {radius_mag * std::cos(convertRadians(true_anomaly))},
        {radius_mag * std::sin(convertRadians(true_anomaly))},
        {0}
    };

    std::vector<std::vector<double>> satellite_position_vector = multiply_two_matrices(transformation_matrix, radius_ijk);

    return satellite_position_vector;

}

std::vector<std::vector<double>> multiply_two_matrices(const std::vector<std::vector<double>> &matrix1, const std::vector<std::vector<double>> &matrix2){
    /**
     * Multiplies two matrices which are modelled as 2 dimensional vectors.
     */
    size_t matrix1_rows = matrix1.size();
    size_t matrix1_cols = matrix1[0].size();

    size_t matrix2_rows = matrix2.size();
    size_t matrix2_cols = matrix2[0].size();

    if (matrix1_cols != matrix2_rows){
        throw std::invalid_argument("Matrix dimensions dont match");
    }

    std::vector<std::vector<double>> result(matrix1_rows, std::vector<double>(matrix2_cols, 0.0));

    for (size_t i {0}; i < matrix1_rows; ++i){
        for (size_t j {0}; j < matrix2_cols; ++j){
            for(size_t k{0}; k < matrix1_cols; ++k){
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    return result;
}

std::vector<std::vector<double>> calculate_ground_position_vector(const double latitude, const double longitude){
    /**
     * Calculates position vector from center of earth to given latitude and longitude.
     */
    double agi = 120.7; // From Astronomers Almanac, recorded on Jan 21, 2019 8:00:16.7646 UTC
    double agi_epoch = ((2019 - 1970) * 31536000) + (20 * 86400) + (8 * 3600) + 16.7646;
    std::chrono::system_clock::time_point tp = std::chrono::system_clock::now();
    std::chrono::system_clock::duration dtn = tp.time_since_epoch();
    double current_epoch = dtn.count() * std::chrono::system_clock::period::num / std::chrono::system_clock::period::den;

    double delta_t = current_epoch - agi_epoch;

    double earth_rotations_per_second = EARTH_ROTATIONS_PER_DAY / 86400;

    double ag = std::fmod(agi + (360 * (earth_rotations_per_second * delta_t)), 360);
    double a = ag + longitude;

    std::vector<std::vector<double>> ground_position_vector = {
        {EARTH_RADIUS * std::cos(convertRadians(latitude)) * std::cos(convertRadians(a)), EARTH_RADIUS * std::cos(convertRadians(latitude)) * std::sin(convertRadians(a)), EARTH_RADIUS * std::sin(convertRadians(latitude))}
    };

    return ground_position_vector;
}

double calculate_latency_distance(const std::vector<std::vector<double>> &satellite_position_vector, 
    const std::vector<std::vector<double>> &ground_position_vector){
    /**
     * Calculates the distance between a ground station and a satellite given the position vectors 
     * to both from the center of the earth.
     */
    std::vector<std::vector<double>> latency_vector = {
        {satellite_position_vector[0][0] - ground_position_vector[0][0], satellite_position_vector[0][1] - ground_position_vector[0][1], satellite_position_vector[0][2] - ground_position_vector[0][2]}
    };

    double latency_distance = std::sqrt((latency_vector[0][0] * latency_vector[0][0]) + (latency_vector[0][1] * latency_vector[0][1]) + (latency_vector[0][2] * latency_vector[0][2]));
    return latency_distance;
}

void display_results(const double true_anomaly, const double latency_distance, const double latency){
    /**
     * Displays final calculated true_anomaly, latency distance, and latency.
     */
    std::cout << std::endl << "---------------------------------------" << std::endl;
    std::cout << "Calculated Orbital Data: " << std::endl << std::endl;
    std::cout << "True Anomaly: " << true_anomaly << "°" << std::endl;
    std::cout << "Satellite-Ground Distance: " << latency_distance << " m" << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "Latency: " << latency << " ms" << std::endl;
    std::cout << "---------------------------------------" << std::endl;
}


int main(){
    std::string satellite_name {};
    std::string TLE_line1 {};
    std::string TLE_line2 {};

    std::cin.clear();
    std::cout << "Enter TLE data for satellite (3 lines, (satellite name on first line and TLE data on next two)): " << std::endl;
    std::getline(std::cin, satellite_name);
    std::getline(std::cin, TLE_line1);
    std::getline(std::cin, TLE_line2);

    SatelliteMetadata metadata = parseTLEline1(TLE_line1);
    OrbitalParams params = parseTLEline2(TLE_line2);

    double latitude {};
    double longitude {};

    std::cin.clear();
    std::cout << "Enter your current latitude in degrees: " << std::endl;
    std::cin >> latitude;
    std::cout << "Enter your current longitude in degrees: " << std::endl;
    std::cin >> longitude;

    char selection {};
    do {
        display_data(satellite_name, metadata, params, longitude, latitude);
        std::cout << "Would you like to calculate latency for the following data? (Y/n) " << std::endl;
        std::cin >> selection;
        selection = toupper(selection);
        if (selection != 'Y' && selection != 'N'){
            std::cout << "Invalid input. Try again." << std::endl;
        }
    } while (selection != 'Y' && selection != 'N');

    if (selection == 'Y'){
         // display_data(satellite_name, metadata, params, latitude, longitude);
         double true_anomaly = calculate_true_anomaly(params.mean_motion,params.mean_anomaly, metadata.epoch_year, metadata.epoch_day, params.eccentricity);
         std::vector<std::vector<double>> satellite_position_vector = calculate_satellite_position_vector(params.inclination, params.right_ascension, params.argument_of_perigee, params.mean_motion, true_anomaly, params.eccentricity);
         std::vector<std::vector<double>> ground_position_vector = calculate_ground_position_vector(latitude, longitude);
         double latency_distance = calculate_latency_distance(satellite_position_vector, ground_position_vector);
         double latency = ((latency_distance * 2) / SPEED_OF_LIGHT) * 1000; // in ms

         display_results(true_anomaly, latency_distance, latency);

    } else {
        std::cout << "Goodbye" << std::endl;
    }
    return 0;
}