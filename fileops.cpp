#include <cctype>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Parameter {
  int paramID;
  double value;
  double variance;
  std::string name;
};

struct Transition {
  int j_prime;
  int ka_prime;
  int kc_prime;
  int v_prime;
  int j_double_prime;
  int ka_double_prime;
  int kc_double_prime;
  int v_double_prime;
  double freq;
  double uncert;
  double weight;
};

struct PredictedTransition {
  double freq;
  double error;
  double logint;
  int freedom;
  double lower_energy;
  int upper_degeneracy;
  int species;
  int j_prime;
  int ka_prime;
  int kc_prime;
  int v_prime;
  int j_double_prime;
  int ka_double_prime;
  int kc_double_prime;
  int v_double_prime;
};

std::vector<Parameter> readParameters(std::string &filename) {
  int id;
  double val;
  double var;
  std::string name;

  std::vector<std::string> v;
  std::vector<Parameter> par;
  std::string line;

  filename += ".par";
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }
  while (std::getline(file, line)) {
    v.push_back(line);
  }
  for (int i = 3; i < 28; i++) {
    std::istringstream ss(v[i]);
    while (ss >> id >> val >> var >> name) {
      Parameter tmp = {id, val, var, name};
      par.push_back(tmp);
    }
  }
  // loop through par and print out the values
  for (auto &p : par) {
    std::cout << p.paramID << " " << p.value << " " << p.variance << " "
              << p.name << std::endl;
  }

  return par;
}

std::vector<Transition> readTransitions(std::string &filename) {
  std::string j_prime;
  std::string ka_prime;
  std::string kc_prime;
  std::string v_prime;
  std::string j_double_prime;
  std::string ka_double_prime;
  std::string kc_double_prime;
  std::string v_double_prime;
  std::string extra;
  int j_prime_int;
  int ka_prime_int;
  int kc_prime_int;
  int v_prime_int;
  int j_double_prime_int;
  int ka_double_prime_int;
  int kc_double_prime_int;
  int v_double_prime_int;
  std::string freq_str;
  std::string uncert_str;
  std::string weight_str;
  double freq;
  double uncert;
  double weight;

  std::string line;
  std::string empty_line = "";

  std::vector<Transition> transitions;
  std::vector<std::string> v;

  filename += ".lin";
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open .lin file: " + filename);
  }
  while (std::getline(file, line)) {
    if (line != empty_line) {
      std::cout << line << std::endl;
      v.push_back(line);
    }
  }
  for (auto &line : v) {
    j_prime = line.substr(0, 3);
    ka_prime = line.substr(3, 3);
    kc_prime = line.substr(6, 3);
    v_prime = line.substr(9, 3);
    j_double_prime = line.substr(12, 3);
    ka_double_prime = line.substr(15, 3);
    kc_double_prime = line.substr(18, 3);
    v_double_prime = line.substr(21, 3);
    extra = line.substr(24, 12);

    int j_prime_int = std::stoi(j_prime);
    int ka_prime_int = std::stoi(ka_prime);
    int kc_prime_int = std::stoi(kc_prime);
    int v_prime_int = std::stoi(v_prime);
    int j_double_prime_int = std::stoi(j_double_prime);
    int ka_double_prime_int = std::stoi(ka_double_prime);
    int kc_double_prime_int = std::stoi(kc_double_prime);
    int v_double_prime_int = std::stoi(v_double_prime);

    double freq = std::stod(line.substr(36, 15));
    double uncert = std::stod(line.substr(52, 12));
    double weight = std::stod(line.substr(65, 10));

    Transition tmp = {j_prime_int,
                      ka_prime_int,
                      kc_prime_int,
                      v_prime_int,
                      j_double_prime_int,
                      ka_double_prime_int,
                      kc_double_prime_int,
                      v_double_prime_int,
                      freq,
                      uncert,
                      weight};
    transitions.push_back(tmp);
  }
  for (auto &t : transitions) {
    std::cout << std::scientific << t.j_prime << " " << t.ka_prime << " "
              << t.kc_prime << " " << t.v_prime << " " << t.j_double_prime
              << " " << t.ka_double_prime << " " << t.kc_double_prime << " "
              << t.v_double_prime << " " << t.freq << " " << t.uncert << " "
              << t.weight << std::endl;
  }
  return transitions;
}

std::vector<PredictedTransition>
readcat(std::string &filename, double lower_bound, double upper_bound) {
  std::vector<PredictedTransition> predictions;
  int convert_qn(std::string qn);

  double freq;
  double error;
  double logint;
  int freedom;
  double lower_energy;
  int upper_degeneracy;
  int species;
  std::string qnfmt;
  std::string j_prime;
  std::string ka_prime;
  std::string kc_prime;
  std::string v_prime;
  std::string extra;
  std::string j_double_prime;
  std::string ka_double_prime;
  std::string kc_double_prime;
  std::string v_double_prime;

  std::string line;
  std::string empty_line = "";

  std::vector<Transition> transitions;
  std::vector<std::string> v;

  filename += ".cat";
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open .cat file: " + filename);
  }
  while (std::getline(file, line)) {
    if (line != empty_line) {
      v.push_back(line);
    }
  }

  for (auto &line : v) {
    freq = std::stod(line.substr(0, 13));
    if (freq > lower_bound && freq < upper_bound) {
      error = std::stod(line.substr(13, 8));
      logint = std::stod(line.substr(21, 8));
      freedom = std::stoi(line.substr(29, 2));
      lower_energy = std::stod(line.substr(31, 10));
      upper_degeneracy = std::stoi(line.substr(41, 3));
      species = std::stoi(line.substr(44, 7));
      qnfmt = line.substr(51, 4);

      j_prime = line.substr(55, 2);
      ka_prime = line.substr(57, 2);
      kc_prime = line.substr(59, 2);
      v_prime = line.substr(61, 2);
      extra = line.substr(63, 4);
      j_double_prime = line.substr(67, 2);
      ka_double_prime = line.substr(69, 2);
      kc_double_prime = line.substr(71, 2);
      v_double_prime = line.substr(73, 2);

      int j_prime_int = convert_qn(j_prime);
      int ka_prime_int = convert_qn(ka_prime);
      int kc_prime_int = convert_qn(kc_prime);
      int v_prime_int = convert_qn(v_prime);
      int j_double_prime_int = convert_qn(j_double_prime);
      int ka_double_prime_int = convert_qn(ka_double_prime);
      int kc_double_prime_int = convert_qn(kc_double_prime);
      int v_double_prime_int = convert_qn(v_double_prime);

      PredictedTransition tmp = {freq,
                                 error,
                                 logint,
                                 freedom,
                                 lower_energy,
                                 upper_degeneracy,
                                 species,
                                 j_prime_int,
                                 ka_prime_int,
                                 kc_prime_int,
                                 v_prime_int,
                                 j_double_prime_int,
                                 ka_double_prime_int,
                                 kc_double_prime_int,
                                 v_double_prime_int};
      predictions.push_back(tmp);
    } else {
      continue;
    }
  }
  /*
  for (auto &t : predictions) {
    std::cout << std::scientific << t.freq << " " << t.error << " " << t.logint
              << " " << t.freedom << " " << t.lower_energy << " "
              << t.upper_degeneracy << " " << t.species << " " << t.j_prime
              << " " << t.ka_prime << " " << t.kc_prime << " " << t.v_prime
              << " " << t.j_double_prime << " " << t.ka_double_prime << " "
              << t.kc_double_prime << " " << t.v_double_prime << std::endl;
  }
  */

  return predictions;
}

void write_dummylin(std::string &filename,
                    std::vector<PredictedTransition> &predictions) {
  filename += "dummy.lin";
  std::ofstream file(filename, std::ofstream::out | std::ofstream::trunc);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }
  for (auto &p : predictions) {
    file << std::setw(3) << p.j_prime << std::setw(3) << p.ka_prime
         << std::setw(3) << p.kc_prime << std::setw(3) << p.v_prime
         << std::setw(3) << p.j_double_prime << std::setw(3) << std::setw(3)
         << p.ka_double_prime << std::setw(3) << p.kc_double_prime
         << std::setw(3) << p.v_double_prime << std::setw(3) << 0
         << std::setw(3) << 0 << std::setw(3) << 0 << std::setw(3) << 0
         << std::setw(16) << std::setprecision(6) << std::fixed << p.freq
         << std::setw(12) << 0.050000 << std::setw(10) << 1.000000 << std::endl;
  }
}

int convert_qn(std::string &qn) {
  char first_letter = qn.at(0);
  int tmp{0};
  if (isalpha(first_letter)) {
    for (auto &c : qn) {
      if (isupper(c)) {
        tmp += 10 * ((int)c - 55);
      } else if (islower(c)) {
        tmp += -10 * ((int)c - 96);
      } else {
        tmp += (int)c - '0';
      }
    }
  } else {
    return std::stoi(qn);
  }
  return tmp;
}

int main() {
  std::string filename = "cyanomethcycloprop_gs";
  std::vector<PredictedTransition> predictions =
      readcat(filename, 250000, 500000);
  write_dummylin(filename, predictions);
}
