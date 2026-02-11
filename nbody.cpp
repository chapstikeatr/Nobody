#include <chrono> 
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

const int DIM = 3;
const double GRAV = 6.67408e-11;
const double SF = 0.01; // softening factor
struct particle {
  double mass;
  vector<double> velocity; // size 3
  vector<double> position; // size 3
  vector<double> force;    // size 3
};

static bool is_integer_string(const string &s) {
  if (s.empty())
    return false;
  size_t i = 0;
  if (s[0] == '+' || s[0] == '-')
    i = 1;
  if (i >= s.size())
    return false;
  for (; i < s.size(); i++) {
    if (!isdigit(static_cast<unsigned char>(s[i])))
      return false;
  }
  return true;
}

static vector<string> parse_tsv(const string &line) {
  vector<string> out;
  string cell;
  istringstream iss(line);
  while (getline(iss, cell, '\t'))
    out.push_back(cell);
  return out;
}

static vector<particle> load_state(string &path) {
  ifstream f(path);
  if (!f)
    throw runtime_error("Could not open input file: " + path);

  string line;
  while (getline(f, line)) {
    bool all_ws = true;
    for (char c : line) {
      if (!isspace(static_cast<unsigned char>(c))) {
        all_ws = false;
        break;
      }
    }
    if (!all_ws)
      break;
  }
  if (line.empty())
    throw runtime_error("No state line found in file.");

  auto cells = parse_tsv(line);
  if (cells.empty())
    throw runtime_error("Bad TSV line.");

  int n = stoi(cells[0]);
  if (n <= 0)
    throw runtime_error("Particle count must be <= 0");

  int DPV = 10; // mass + 3 pos + 3 vel + 3 force
  if ((int)cells.size() < 1 + n * DPV)
    throw runtime_error(
        "TSV line does not have enough columns for n particles.");

  vector<particle> ps;
  ps.reserve((size_t)n);

  size_t idx = 1;
  for (int i = 0; i < n; i++) {
    particle p;
    p.velocity.assign(3, 0.0);
    p.position.assign(3, 0.0);
    p.force.assign(3, 0.0);

    p.mass = stod(cells[idx++]);

    p.position[0] = stod(cells[idx++]);
    p.position[1] = stod(cells[idx++]);
    p.position[2] = stod(cells[idx++]);

    p.velocity[0] = stod(cells[idx++]);
    p.velocity[1] = stod(cells[idx++]);
    p.velocity[2] = stod(cells[idx++]);

    p.force[0] = stod(cells[idx++]);
    p.force[1] = stod(cells[idx++]);
    p.force[2] = stod(cells[idx++]);

    ps.push_back(p);
  }

  return ps;
}

static vector<particle> random_init(int n) {
  if (n <= 0)
    throw runtime_error("n must be > 0");

  srand((unsigned)time(NULL));

  vector<particle> ps;
  ps.reserve((size_t)n);

  for (int i = 0; i < n; i++) {
    particle p;
    p.velocity.assign(3, 0.0);
    p.position.assign(3, 0.0);
    p.force.assign(3, 0.0);

    // simple demo ranges; tune as desired
    p.mass = 1e20 + (double)rand() / RAND_MAX * (1e30 - 1e20);
    p.position[0] = -1.5e11 + (double)rand() / RAND_MAX * (3.0e11);
    p.position[1] = -1.5e11 + (double)rand() / RAND_MAX * (3.0e11);
    p.position[2] = -1.5e11 + (double)rand() / RAND_MAX * (3.0e11);
    p.velocity[0] = -3e4 + (double)rand() / RAND_MAX * (6.0e4);
    p.velocity[1] = -3e4 + (double)rand() / RAND_MAX * (6.0e4);
    p.velocity[2] = -3e4 + (double)rand() / RAND_MAX * (6.0e4);

    ps.push_back(p);
  }

  return ps;
}

void force_calculation(vector<particle> &particles) {
  int np = particles.size();
  double eps2 = SF * SF;

  // reset forces
  for (int i = 0; i < np; i++) {
    particles[i].force[0] = particles[i].force[1] = particles[i].force[2] = 0.0;
  }

  for (int i = 0; i < np; i++) {
    for (int j = i + 1; j < np; j++) {
      double dx = particles[j].position[0] - particles[i].position[0];
      double dy = particles[j].position[1] - particles[i].position[1];
      double dz = particles[j].position[2] - particles[i].position[2];

      double r2 = dx * dx + dy * dy + dz * dz + eps2;
      double inv_r = 1.0 / sqrt(r2);
      double inv_r3 = inv_r * inv_r * inv_r;

      double s = GRAV * particles[i].mass * particles[j].mass * inv_r3;

      double fx = s * dx;
      double fy = s * dy;
      double fz = s * dz;

      particles[i].force[0] += fx;
      particles[i].force[1] += fy;
      particles[i].force[2] += fz;

      particles[j].force[0] -= fx;
      particles[j].force[1] -= fy;
      particles[j].force[2] -= fz;
    }
  }
}

void update_pos(vector<particle> &particles, double dt) {
  for (auto &p : particles) {
    double ax = p.force[0] / p.mass;
    double ay = p.force[1] / p.mass;
    double az = p.force[2] / p.mass;

    p.velocity[0] = p.velocity[0] + ax * dt;
    p.velocity[1] = p.velocity[1] + ay * dt;
    p.velocity[2] = p.velocity[2] + az * dt;

    p.position[0] = p.position[0] + p.velocity[0] * dt;
    p.position[1] = p.position[1] + p.velocity[1] * dt;
    p.position[2] = p.position[2] + p.velocity[2] * dt;
  }
}

// Dumps one state line
static void dump_state(ostream &out, const vector<particle> &ps) {
  out << ps.size();
  out << setprecision(17);
  for (const auto &p : ps) {
    out << '\t' << p.mass << '\t' << p.position[0] << '\t' << p.position[1]
        << '\t' << p.position[2] << '\t' << p.velocity[0] << '\t'
        << p.velocity[1] << '\t' << p.velocity[2] << '\t' << p.force[0] << '\t'
        << p.force[1] << '\t' << p.force[2];
  }
  out << '\n';
}

static void usage(const char *prog) {
  cerr << "Usage:\n"
       << "  " << prog
       << " <N|input.tsv> <dt> <steps> <dump_every> [output.tsv]\n\n"
       << "Examples:\n"
       << "  " << prog << " solar.tsv 200 5000000 1000 solar_out.tsv\n"
       << "  " << prog << " 100 1 10000 10 out.tsv\n";
}

int main(int argc, char *argv[]) {
  try {
    if (argc < 5) {
      usage(argv[0]);
      return 1;
    }

    string model = argv[1];
    double dt = stod(argv[2]);
    int steps = stoll(argv[3]);
    int dump_every = stoll(argv[4]);

    if (dt <= 0)
      throw runtime_error("dt must be > 0");
    if (steps <= 0)
      throw runtime_error("steps must be > 0");
    if (dump_every <= 0)
      throw runtime_error("dump_every must be > 0");

    vector<particle> particles;
    if (is_integer_string(model)) {
      particles = random_init(stoi(model));
    } else {
      particles = load_state(model);
    }

    ofstream fout;
    ostream *out = &cout;
    if (argc >= 6) {
      fout.open(argv[5]);
      if (!fout)
        throw runtime_error("Could not open output file for writing.");
      out = &fout;
    }
    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();

    for (int step = 0; step <= steps; step++) {
      force_calculation(particles);

      if (step % dump_every == 0) {
        dump_state(*out, particles);
      }

      update_pos(particles, dt);
    }

    auto t1 = clock::now();
    std::chrono::duration<double> elapsed = t1 - t0;

    std::cerr << "TIMING seconds=" << elapsed.count()
              << " N=" << particles.size();
    return 0;
  } catch (const exception &e) {
    cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}
