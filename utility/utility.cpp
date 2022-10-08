#include "utility.hpp"

namespace Utility
{
    // https://stackoverflow.com/a/4609795
    template <typename T>
    int sgn(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

    int mod(int a, int b) { return ((a % b) + b) % b; }

    void getParticlesFromCSV(std::string file, std::vector<Particle> &particles)
    {
        rapidcsv::Document doc(file);
        std::vector<double> row;
        for (size_t i = 0; i < doc.GetRowCount(); i++) {
            row = doc.GetRow<double>(i);
            int ID = doc.GetCell<int>(0, i);
            particles.push_back(
                Particle(ID, row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10]));
        }
    }

    void getParticlesFromTuple(
        std::vector<std::tuple<int, double, double, double, double, double, double, double, double, double, double>>
            &tuples,
        std::vector<Particle> &particles)
    {
        for (std::tuple<int, double, double, double, double, double, double, double, double, double, double> &t :
             tuples) {
            particles.push_back(Particle(std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t), std::get<4>(t),
                                         std::get<5>(t), std::get<6>(t), std::get<7>(t), std::get<8>(t), std::get<9>(t),
                                         std::get<10>(t)));
        }
    }

    void writeStepToCSV(std::string file, std::vector<Particle> &particles)
    {
        std::ofstream csvFile;
        csvFile.open(file);
        csvFile << "ID, pX, pY, pZ, vX, vY, vZ, aX, aY, aZ, m\n";
        // std::cout << particles.size();
        for (size_t i = 0; i < particles.size(); i++) {
            if (particles[i].isDummy) continue;
            csvFile << particles[i].ID << ", " << particles[i].pX << ", " << particles[i].pY << ", " << particles[i].pZ
                    << ", " << particles[i].vX << ", " << particles[i].vY << ", " << particles[i].vZ << ", "
                    << particles[i].aX << ", " << particles[i].aY << ", " << particles[i].aZ << ", "
                    << particles[i].mass << "\n";
        }
        csvFile.close();
    }

    void writeStepToCSVWithForces(std::string file, std::vector<Particle> &particles)
    {
        std::ofstream csvFile;
        csvFile.open(file);
        csvFile << "ID, pX, pY, pZ, vX, vY, vZ, aX, aY, aZ, m, fX, fY, fZ\n";
        // std::cout << particles.size();
        for (size_t i = 0; i < particles.size(); i++) {
            if (particles[i].isDummy) continue;
            csvFile << particles[i].ID << ", " << particles[i].pX << ", " << particles[i].pY << ", " << particles[i].pZ
                    << ", " << particles[i].vX << ", " << particles[i].vY << ", " << particles[i].vZ << ", "
                    << particles[i].aX << ", " << particles[i].aY << ", " << particles[i].aZ << ", "
                    << particles[i].mass << particles[i].fX << ", " << particles[i].fY << ", " << particles[i].fZ
                    << "\n";
        }
        csvFile.close();
    }

    // https://stackoverflow.com/a/55422097
    int BinomialCoefficient(const int n, const int k)
    {
        std::vector<int> aSolutions(k);
        aSolutions[0] = n - k + 1;

        for (int i = 1; i < k; ++i) {
            aSolutions[i] = aSolutions[i - 1] * (n - k + 1 + i) / (i + 1);
        }

        return aSolutions[k - 1];
    }

    // https://stackoverflow.com/a/69552155
    std::string get_file_contents(const char *filename)
    {
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        if (!in) {
            std::cerr << "could not open " << filename << std::endl;
            exit(1);
        }
        std::ostringstream contents;
        contents << in.rdbuf();
        return contents.str();
    }

    std::vector<int> getDecomposition(int worldSize,
                                      const std::vector<std::pair<int, std::vector<int>>> &decompositions)
    {
        std::vector<int> result;
        for (std::pair<int, std::vector<int>> e : decompositions) {
            if (e.first == worldSize) {
                result = e.second;
                break;
            }
        }
        return result;
    }
    /*
        std::tuple<std::vector<std::pair<int, std::vector<int>>>, std::vector<std::pair<int, std::vector<int>>>>
        jsonToDecompositionArray(const char *filename, std::string type)
        {
            // std::string json = get_file_contents(filename);
            const char *json = get_file_contents(filename).c_str();

            rapidjson::Document document;
            document.Parse(json);

            // rapidjson::StringBuffer buffer;
            // rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
            // document.Accept(writer);
            // std::cout << buffer.GetString() << std::endl;

            std::tuple<std::vector<std::pair<int, std::vector<int>>>, std::vector<std::pair<int, std::vector<int>>>>
                decompositions;

            const rapidjson::Value &optimal = document["optimal"];
            const rapidjson::Value &naive = document["naive"];

            decompositions = std::tuple(std::vector<std::pair<int, std::vector<int>>>(),
                                        std::vector<std::pair<int, std::vector<int>>>());

            for (rapidjson::Value::ConstMemberIterator it = optimal.MemberBegin(); it != optimal.MemberEnd(); ++it) {
                std::pair<int, std::vector<int>> decomposition;
                decomposition.second = std::vector<int>();
                decomposition.first = std::stoi(it->name.GetString());
                for (rapidjson::SizeType i = 0; i < it->value.Size(); i++) {
                    decomposition.second.push_back(it->value[i].GetInt());
                }
                std::get<0>(decompositions).push_back(decomposition);
            }

            for (rapidjson::Value::ConstMemberIterator it = naive.MemberBegin(); it != naive.MemberEnd(); ++it) {
                std::pair<int, std::vector<int>> decomposition;
                decomposition.second = std::vector<int>();
                decomposition.first = std::stoi(it->name.GetString());
                for (rapidjson::SizeType i = 0; i < it->value.Size(); i++) {
                    decomposition.second.push_back(it->value[i].GetInt());
                }
                std::get<1>(decompositions).push_back(decomposition);
            }

            return decompositions;
        }
    */

}  // namespace Utility