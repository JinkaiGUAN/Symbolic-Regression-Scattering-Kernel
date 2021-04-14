#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace std;

int main()
{
    string line, fName;
    int id, n, i, t, tTime, nAtoms, yStep, currentTimeStep, typ, bo;

    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, r;
    double xLo, xHi, yLo, yHi, zLo, zHi, Lz;

    double TAU1Left, TAU2Left, TAU3Left, TAU1Right, TAU2Right, TAU3Right, keRight;
    double pistonLeftX, pistonRightX, barrierLeftX, barrierRightX;
    int pistonLeftCount, pistonRightCount, barrierLeftCount, barrierRightCount;
    int nMethaneLeft, nMethaneRight, nMethaneMiddle;
    double heightLeft, heightRight, volumeLeft, volumeRight, rhoLeft, rhoRight, pLeft, pRight, tempRight;
    double TAU1Middle, TAU2Middle, TAU3Middle, keMiddle, pMiddle, tempMiddle;
    double xMin, xMax;
    double kB = 1.38064852e-23;

    // Parameters which should be consistent with LAMMPS file.
    int tSkip = 20;
    double deltaT = 2.0;

    double mi = 2.66389e-26; // mass of one molecule
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;

    // ***********************
    double Tm = 600;                    // initial temperature of methane
    double Tw = 423;                    // temperature of wall (Kelvin)
    double vM = sqrt(2 * kB * Tw / mi); // most probable speed
    double vF = sqrt(kB * Tm / mi);     // reference velocity factor
    cout << "Most probable speed is: " << vM << " m/s" << endl;
    cout << "Reference velocity factor is: " << vF << endl;

    //**** THINGS TO CHANGE *****

    double d = 4.5; // should be consistent with combine.cpp

    double rCut = 20;
    double spare = rCut * 0.5;
    double midpoint = rCut + spare + d;


    int nTimeSteps = 401; // CHANGE, use command line: grep -o 'TIMESTEP' dump_meas.lammpstrj | wc -l

    double binWidth = 10; // binwidth of velocity
    int maxVout = vM * 4; // max output velocity
    int nBins = ceil(maxVout / binWidth);
    int tBins = ceil(nTimeSteps / binWidth);

    // **************************

    ifstream data("dump_meas.lammpstrj", ios::in);

    ofstream file3("3_meas_velocities_top.txt", ios::out);
    ofstream file4("3_meas_velocities_bottom.txt", ios::out);

    vector<int> crossedStart;

    vector<double> vxStart;
    vector<double> vyStart;
    vector<double> vzStart;

    vector<double> vxEnd;
    vector<double> vyEnd;
    vector<double> vzEnd;

    vector<double> keStart;
    vector<double> keEnd;

    vector<int> tStart;
    vector<int> tEnd;

    int nFinishes = 0;
    int nStarts = 0;

    for (t = 0; t < nTimeSteps; t++)
    {

        // applicable to all the post-process of LAMMPS files in the future.
        for (n = 1; n < 10; n++)
        {
            if (n == 4)
            {
                data >> nAtoms;
                //                 cout << "nAtoms = " << nAtoms << endl;
            }

            if (n == 2)
            {
                data >> currentTimeStep;

                cout << "currentTimeStep = " << currentTimeStep
                     << "; t = " << t << " [ "
                     << 100 * float(t + 1) / float(nTimeSteps)
                     << "% ]" << endl;
            }

            if (n == 6)
            {
                data >> xLo >> xHi;
            }

            if (n == 7)
            {
                data >> yLo >> yHi;
            }

            if (n == 8)
            {
                data >> zLo >> zHi;
            }

            getline(data, line);
        }

        // set fields
        if (t == 0)
        {
            crossedStart.resize(nAtoms, 0);
            vxStart.resize(nAtoms, 0.0);
            vyStart.resize(nAtoms, 0.0);
            vzStart.resize(nAtoms, 0.0);

            vxEnd.resize(nAtoms, 0.0);
            vyEnd.resize(nAtoms, 0.0);
            vzEnd.resize(nAtoms, 0.0);

            keStart.resize(nAtoms, 0.0);
            keEnd.resize(nAtoms, 0.0);

            tStart.resize(nAtoms, 0);
            tEnd.resize(nAtoms, 0);
        }

        for (n = 0; n < nAtoms; n++)
        {
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> KE >> PE;

            if (typ == 1) // methane
            {
                if (t == 0)
                {
                    if (vz > 0)
                    {
                        crossedStart[n] = 1; // going upwards

                        vxStart[n] = vx * refVelocity;
                        vyStart[n] = vy * refVelocity;
                        vzStart[n] = vz * refVelocity;
                        keStart[n] = KE * refVelocity * refVelocity; 
                    }

                    if (vz < 0)
                    {
                        crossedStart[n] = -1; // going downwards

                        vxStart[n] = vx * refVelocity;
                        vyStart[n] = vy * refVelocity;
                        vzStart[n] = vz * refVelocity;
                        keStart[n] = KE * refVelocity * refVelocity;
                    }
                }
                else
                {
                    if ((z < (midpoint - spare - d)) && (crossedStart[n] == -1)) // downwards
                    {
                        nStarts++;
                        crossedStart[n] = -2; // started
                        tStart[n] = t;
                    }
                    else if ((z > (midpoint - spare - d)) && (crossedStart[n] == -2)) // finished
                    {
                        crossedStart[n] = -3;
                        vxEnd[n] = vx * refVelocity;
                        vyEnd[n] = vy * refVelocity;
                        vzEnd[n] = vz * refVelocity;
                        nFinishes++;
                        tEnd[n] = t;
                    }

                    if ((z > (midpoint + spare + d)) && (crossedStart[n] == 1)) // upwards
                    {
                        nStarts++;
                        crossedStart[n] = 2; //started
                        tStart[n] = t;
                    }
                    else if ((z < (midpoint + spare + d)) && (crossedStart[n] == 2)) // finished
                    {
                        crossedStart[n] = 3;
                        vxEnd[n] = vx * refVelocity;
                        vyEnd[n] = vy * refVelocity;
                        vzEnd[n] = vz * refVelocity;
                        nFinishes++;
                        tEnd[n] = t;
                    }
                }
            }
        }

        getline(data, line);

        cout << endl
             << "number of starts = " << nStarts << ", number of finishes = " << nFinishes
             << ", as a percentage of total atoms: " << 100 * nFinishes / nAtoms << "%"
             << endl;
        cout << endl;
    }

    
    // dump the incident and reflected velocity in raw data.
    vector<double> vN, vNi;
    vector<double> vT, vTi;
    vector<double> vTx, vTxi;
    vector<double> vTy, vTyi;
    double vMagT, vMagTi;
    vector<int> residenceT;

    for (i = 0; i < nAtoms; i++)
    {
        if (crossedStart[i] == 3)
        {
            file3 << vxStart[i] << " " << vyStart[i] << " " << vzStart[i] << " "
                  << vxEnd[i] << " " << vyEnd[i] << " " << vzEnd[i] << endl;

            // normal
            vN.push_back(abs(vzEnd[i]));
            vNi.push_back(abs(vzStart[i]));

            // combined tangents
            vMagT = vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i];
            vMagTi = vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i];

            if (vMagT > 0)
            {
                vT.push_back(sqrt(vMagT));
            }
            else
            {
                vT.push_back(0.0);
            }

            if (vMagTi > 0)
            {
                vTi.push_back(sqrt(vMagTi));
            }
            else
            {
                vTi.push_back(0.0);
            }

            vTx.push_back(vxEnd[i]);
            vTxi.push_back(vxStart[i]);

            vTy.push_back(vyEnd[i]);
            vTyi.push_back(vyStart[i]);
            residenceT.push_back(tEnd[i] - tStart[i]);
        }

        if (crossedStart[i] == -3)
        {
            file4 << vxStart[i] << " " << vyStart[i] << " " << vzStart[i] << " "
                  << vxEnd[i] << " " << vyEnd[i] << " " << vzEnd[i] << endl;

            // normal
            vN.push_back(abs(vzEnd[i]));
            vNi.push_back(abs(vzStart[i]));

            // combined tangents
            vMagT = vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i];
            vMagTi = vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i];

            if (vMagT > 0)
            {
                vT.push_back(sqrt(vMagT));
            }
            else
            {
                vT.push_back(0.0);
            }

            if (vMagTi > 0)
            {
                vTi.push_back(sqrt(vMagTi));
            }
            else
            {
                vTi.push_back(0.0);
            }

            vTx.push_back(vxEnd[i]);
            vTxi.push_back(vxStart[i]);

            vTy.push_back(vyEnd[i]);
            vTyi.push_back(vyStart[i]);
            residenceT.push_back(tEnd[i] - tStart[i]);
        }

        // residenceT.push_back(tEnd[i] - tStart[i]);
    }

    int nPts = vN.size();

    double sum_vTx = 0.0, sum_vTxi = 0.0, ave_vTx, ave_vTxi;
    for (int i = 0; i < vTxi.size(); ++i)
    {
        sum_vTxi += vTxi[i];
    }
    for (int i = 0; i < vTx.size(); ++i)
    {
        sum_vTx += vTx[i];
    }
    ave_vTxi = sum_vTxi / vTxi.size();
    ave_vTx = sum_vTx / vTx.size();

    cout << "The average incident velocity of tangential X components: " << ave_vTxi << endl;
    cout << "The average reflected velocity of tangential X components: " << ave_vTx << endl;

    // Residence time histogram
    {
        vector<double> distribution(tBins, 0.0);
        vector<double> bins(tBins, 0.0);

        for (i = 0; i < tBins; i++)
        {
            bins[i] = binWidth * 0.5 + binWidth * i;
        }

        for (i = 0; i < nPts; i++)
        {
            bo = floor(residenceT[i] / binWidth);

            if (bo >= tBins)
            {
                bo = tBins - 1;
            }

            distribution[bo] += 1.0;
        }

        vector<double> probability(tBins, 0.0);

        // scale with area
        double areaUnderGraph = 0.0;

        for (i = 0; i < tBins - 1; i++)
        {
            areaUnderGraph += (bins[i + 1] - bins[i]) * (distribution[i] + distribution[i + 1]) * 0.5;
        }

        for (i = 0; i < tBins; i++)
        {
            probability[i] = distribution[i] / areaUnderGraph;
        }

        ofstream file("residence_time_distribution.txt");

        for (i = 0; i < tBins; i++)
        {
            file << bins[i] << " " << distribution[i] << " " << probability[i] << endl;
        }
    }

    // ************************ incident normal component *******************************
    {
        // now let's output some histogram
        double vMax = 0.0;
        vector<double> distribution(nBins, 0.0);
        vector<double> bins(nBins, 0.0);

        for (i = 0; i < nBins; i++)
        {
            bins[i] = binWidth * 0.5 + binWidth * i;
        }

        for (i = 0; i < nPts; i++)
        {
            if (vNi[i] < maxVout)
            {
                bo = floor(vNi[i] / binWidth);

                if (bo >= nBins)
                {
                    bo = nBins - 1;
                }

                distribution[bo] += 1.0;
            }
            else
            {
                cout << "Warning: value out of bounds for distribution, set distribution larger "
                     << "Value of measure velocity = " << vNi[i]
                     << endl;
            }

            if (vNi[i] > vMax)
            {
                vMax = vNi[i];
            }
        }

        cout << "Maximum vNi value = " << vMax << endl;

        vector<double> probability(nBins, 0.0);

        // scale with area
        double areaUnderGraph = 0.0;

        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (bins[i + 1] - bins[i]) * (distribution[i] + distribution[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probability[i] = distribution[i] / areaUnderGraph;
        }

        ofstream file("incident_velocity_normal.txt");

        for (i = 0; i < nBins; i++)
        {
            file << bins[i] << " " << distribution[i] << " " << probability[i] << endl;
        }
    }
    // **********************************************************************************

    // ************************ reflected normal component ******************************
    {
        // now let's output some histogram
        double vMax = 0.0;
        vector<double> distribution(nBins, 0.0);
        vector<double> bins(nBins, 0.0);

        for (i = 0; i < nBins; i++)
        {
            bins[i] = binWidth * 0.5 + binWidth * i;
        }

        for (i = 0; i < nPts; i++)
        {

            if (vN[i] < maxVout)
            {
                bo = floor(vN[i] / binWidth);

                if (bo >= nBins)
                {
                    bo = nBins - 1;
                }

                distribution[bo] += 1.0;
            }
            else
            {
                cout << "Warning: value out of bounds of distribution, set distribution larger "
                     << "Value of measured velocity = " << vN[i]
                     << endl;
            }

            if (vN[i] > vMax)
            {
                vMax = vN[i];
            }
        }

        cout << "Maximum vN value = " << vMax << endl;

        vector<double> probability(nBins, 0.0);

        // scale with area
        double areaUnderGraph = 0.0;

        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (bins[i + 1] - bins[i]) * (distribution[i] + distribution[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probability[i] = distribution[i] / areaUnderGraph;
        }

        ofstream file("reflected_velocity_normal.txt");

        for (i = 0; i < nBins; i++)
        {
            file << bins[i] << " " << distribution[i] << " " << probability[i] << endl;
        }
    }
    // **********************************************************************************

    //************************* incident tangential combined ****************************
    {
        // now let's output some histogram
        double vMax = 0.0;
        vector<double> distribution(nBins, 0.0);
        vector<double> bins(nBins, 0.0);

        for (i = 0; i < nBins; i++)
        {
            bins[i] = binWidth * 0.5 + binWidth * i;
        }

        for (i = 0; i < nPts; i++)
        {
            if (vTi[i] < maxVout)
            {
                bo = floor(vTi[i] / binWidth);

                if (bo >= nBins)
                {
                    bo = nBins - 1;
                }

                distribution[bo] += 1.0;
            }
            else
            {
                cout << "Warning: value out of bounds of distribution, set distribution larger "
                     << "Value of measured velocity = " << vTi[i]
                     << endl;
            }

            if (vTi[i] > vMax)
            {
                vMax = vTi[i];
            }
        }

        cout << "Maximum vTi value = " << vMax << endl;

        vector<double> probability(nBins, 0.0);

        // scale with area
        double areaUnderGraph = 0.0;

        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (bins[i + 1] - bins[i]) * (distribution[i] + distribution[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probability[i] = distribution[i] / areaUnderGraph;
        }

        ofstream file("incident_velocity_tangential.txt");

        for (i = 0; i < nBins; i++)
        {
            file << bins[i] << " " << distribution[i] << " " << probability[i] << endl;
        }
    }
    //***********************************************************************************

    // ************************ reflected tangential combined ***************************
    {
        // now let's output some histogram
        double vMax = 0.0;
        vector<double> distribution(nBins, 0.0);
        vector<double> bins(nBins, 0.0);

        for (i = 0; i < nBins; i++)
        {
            bins[i] = binWidth * 0.5 + binWidth * i;
        }

        for (i = 0; i < nPts; i++)
        {
            if (vT[i] < maxVout)
            {
                bo = floor(vT[i] / binWidth);

                if (bo >= nBins)
                {
                    bo = nBins - 1;
                }

                distribution[bo] += 1.0;
            }
            else
            {
                cout << "Warning: value out of bounds of distribution, set distribution larger "
                     << "Value of measured velocity = " << vT[i]
                     << endl;
            }

            if (vT[i] > vMax)
            {
                vMax = vT[i];
            }
        }

        cout << "Maximum vT value = " << vMax << endl;

        vector<double> probability(nBins, 0.0);

        // scale with area
        double areaUnderGraph = 0.0;

        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (bins[i + 1] - bins[i]) * (distribution[i] + distribution[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probability[i] = distribution[i] / areaUnderGraph;
        }

        ofstream file("reflected_velocity_tangential.txt");

        for (i = 0; i < nBins; i++)
        {
            file << bins[i] << " " << distribution[i] << " " << probability[i] << endl;
        }
    }
    //***********************************************************************************

    //************************* incident tangential X ***********************************
    {
        // now let's output some histogram
        double vMax = 0.0;
        double vMin = 0.0;

        vector<double> distributionA(nBins, 0.0);
        vector<double> distributionB(nBins, 0.0);
        vector<double> binsA(nBins, 0.0);
        vector<double> binsB(nBins, 0.0);

        for (i = 0; i < nBins; i++)
        {
            binsA[i] = binWidth * 0.5 + binWidth * i;
            binsB[i] = -binsA[i];
        }

        for (i = 0; i < nPts; i++)
        {
            if (abs(vTxi[i]) < maxVout)
            {
                if (vTxi[i] >= 0)
                {
                    bo = floor(vTxi[i] / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionA[bo] += 1.0;
                }
                else
                {
                    bo = floor(abs(vTxi[i]) / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionB[bo] += 1.0;
                }
            }
            else
            {
                cout << "Warning: value out of bounds of distribution, set distribution larger "
                     << "Value of measured velocity = " << vTxi[i]
                     << endl;
            }

            if (vTxi[i] > vMax)
            {
                vMax = vTxi[i];
            }

            if (vTxi[i] < vMin)
            {
                vMin = vTxi[i];
            }
        }

        cout << "Maximum vTxi value = " << vMax << endl;
        cout << "Minimum vTxi value = " << vMin << endl;

        vector<double> probabilityA(nBins, 0.0);
        vector<double> probabilityB(nBins, 0.0);

        // scale with area
        double areaUnderGraph = 0.0;

        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (binsA[i + 1] - binsA[i]) * (distributionA[i] + distributionA[i + 1]) * 0.5;
            areaUnderGraph += abs(binsB[i + 1] - binsB[i]) * (distributionB[i] + distributionB[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probabilityA[i] = distributionA[i] / areaUnderGraph;
            probabilityB[i] = distributionB[i] / areaUnderGraph;
        }

        ofstream file("incident_velocity_tangential_x.txt");

        for (i = 0; i < nBins; i++)
        {
            file << binsB[nBins - i - 1] << " " << distributionB[nBins - i - 1] << " " << probabilityB[nBins - i - 1] << endl;
        }

        for (i = 0; i < nBins; i++)
        {
            file << binsA[i] << " " << distributionA[i] << " " << probabilityA[i] << endl;
        }
    }
    //***********************************************************************************

    //************************* reflected tangential X **********************************
    {
        // now let's output some histogram
        double vMax = 0.0;
        double vMin = 0.0;

        vector<double> distributionA(nBins, 0.0);
        vector<double> distributionB(nBins, 0.0);
        vector<double> binsA(nBins, 0.0);
        vector<double> binsB(nBins, 0.0);

        for (i = 0; i < nBins; i++)
        {
            binsA[i] = binWidth * 0.5 + binWidth * i;
            binsB[i] = -binsA[i];
        }

        for (i = 0; i < nPts; i++)
        {
            if (abs(vTx[i]) < maxVout)
            {
                if (vTx[i] >= 0)
                {
                    bo = floor(vTx[i] / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionA[bo] += 1.0;
                }
                else
                {
                    bo = floor(abs(vTx[i]) / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionB[bo] += 1.0;
                }
            }
            else
            {
                cout << "Warning: value out of bounds of distribution, set distribution larger "
                     << "Value of measured velocity = " << vTx[i]
                     << endl;
            }

            if (vTx[i] > vMax)
            {
                vMax = vTx[i];
            }

            if (vTx[i] < vMin)
            {
                vMin = vTx[i];
            }
        }

        cout << "Maximum vTx value = " << vMax << endl;
        cout << "Minimum vTx value = " << vMin << endl;

        vector<double> probabilityA(nBins, 0.0);
        vector<double> probabilityB(nBins, 0.0);

        // scale with area
        double areaUnderGraph = 0.0;

        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (binsA[i + 1] - binsA[i]) * (distributionA[i] + distributionA[i + 1]) * 0.5;
            areaUnderGraph += abs(binsB[i + 1] - binsB[i]) * (distributionB[i] + distributionB[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probabilityA[i] = distributionA[i] / areaUnderGraph;
            probabilityB[i] = distributionB[i] / areaUnderGraph;
        }

        ofstream file("reflected_velocity_tangential_x.txt");

        for (i = 0; i < nBins; i++)
        {
            file << binsB[nBins - i - 1] << " " << distributionB[nBins - i - 1] << " " << probabilityB[nBins - i - 1] << endl;
        }

        for (i = 0; i < nBins; i++)
        {
            file << binsA[i] << " " << distributionA[i] << " " << probabilityA[i] << endl;
        }
    }
    //***********************************************************************************

    //************************* incident tangential Y ***********************************
    {
        // now let's output some histogram
        double vMax = 0.0;
        double vMin = 0.0;

        vector<double> distributionA(nBins, 0.0);
        vector<double> distributionB(nBins, 0.0);
        vector<double> binsA(nBins, 0.0);
        vector<double> binsB(nBins, 0.0);

        for (i = 0; i < nBins; i++)
        {
            binsA[i] = binWidth * 0.5 + binWidth * i;
            binsB[i] = -binsA[i];
        }

        for (i = 0; i < nPts; i++)
        {
            if (abs(vTyi[i]) < maxVout)
            {
                if (vTyi[i] >= 0)
                {
                    bo = floor(vTyi[i] / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionA[bo] += 1.0;
                }
                else
                {
                    bo = floor(abs(vTyi[i]) / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionB[bo] += 1.0;
                }
            }
            else
            {
                cout << "Warning: value out of bounds of distribution, set distribution larger "
                     << "Value of measured velocity = " << vTyi[i]
                     << endl;
            }

            if (vTyi[i] > vMax)
            {
                vMax = vTyi[i];
            }

            if (vTyi[i] < vMin)
            {
                vMin = vTyi[i];
            }
        }

        cout << "Maximum vTyi value = " << vMax << endl;
        cout << "Minimum vTyi value = " << vMin << endl;

        vector<double> probabilityA(nBins, 0.0);
        vector<double> probabilityB(nBins, 0.0);

        // scale with area
        double areaUnderGraph = 0.0;

        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (binsA[i + 1] - binsA[i]) * (distributionA[i] + distributionA[i + 1]) * 0.5;
            areaUnderGraph += abs(binsB[i + 1] - binsB[i]) * (distributionB[i] + distributionB[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probabilityA[i] = distributionA[i] / areaUnderGraph;
            probabilityB[i] = distributionB[i] / areaUnderGraph;
        }

        ofstream file("incident_velocity_tangential_y.txt");

        for (i = 0; i < nBins; i++)
        {
            file << binsB[nBins - i - 1] << " " << distributionB[nBins - i - 1] << " " << probabilityB[nBins - i - 1] << endl;
        }

        for (i = 0; i < nBins; i++)
        {
            file << binsA[i] << " " << distributionA[i] << " " << probabilityA[i] << endl;
        }
    }
    //***********************************************************************************

    //************************* reflected tangential Y **********************************
    {
        // now let's output some histogram
        double vMax = 0.0;
        double vMin = 0.0;

        vector<double> distributionA(nBins, 0.0);
        vector<double> distributionB(nBins, 0.0);
        vector<double> binsA(nBins, 0.0);
        vector<double> binsB(nBins, 0.0);

        for (i = 0; i < nBins; i++)
        {
            binsA[i] = binWidth * 0.5 + binWidth * i;
            binsB[i] = -binsA[i];
        }

        for (i = 0; i < nPts; i++)
        {
            if (abs(vTy[i]) < maxVout)
            {
                if (vTy[i] >= 0)
                {
                    bo = floor(vTy[i] / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionA[bo] += 1.0;
                }
                else
                {
                    bo = floor(abs(vTy[i]) / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionB[bo] += 1.0;
                }
            }
            else
            {
                cout << "Warning: value out of bounds of distribution, set distribution larger "
                     << "Value of measured velocity = " << vTy[i]
                     << endl;
            }

            if (vTy[i] > vMax)
            {
                vMax = vTy[i];
            }

            if (vTy[i] < vMin)
            {
                vMin = vTy[i];
            }
        }

        cout << "Maximum vTy value = " << vMax << endl;
        cout << "Minimum vTy value = " << vMin << endl;

        vector<double> probabilityA(nBins, 0.0);
        vector<double> probabilityB(nBins, 0.0);

        // scale with area
        double areaUnderGraph = 0.0;

        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (binsA[i + 1] - binsA[i]) * (distributionA[i] + distributionA[i + 1]) * 0.5;
            areaUnderGraph += abs(binsB[i + 1] - binsB[i]) * (distributionB[i] + distributionB[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probabilityA[i] = distributionA[i] / areaUnderGraph;
            probabilityB[i] = distributionB[i] / areaUnderGraph;
        }

        ofstream file("reflected_velocity_tangential_y.txt");

        for (i = 0; i < nBins; i++)
        {
            file << binsB[nBins - i - 1] << " " << distributionB[nBins - i - 1] << " " << probabilityB[nBins - i - 1] << endl;
        }

        for (i = 0; i < nBins; i++)
        {
            file << binsA[i] << " " << distributionA[i] << " " << probabilityA[i] << endl;
        }
    }

    return 0;
}
