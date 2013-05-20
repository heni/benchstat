#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace std;

void apply_permutation(vector<int>& items, const vector<int>& permutation) {
    vector<int> tmp(items.size());
    for (size_t i = 0; i < permutation.size(); ++i) {
        tmp[i] = items[permutation[i]];
    }
    items.swap(tmp);
}


class BaseBenchmark {
protected:
    vector<int> InternalItems;
    vector<int> Permutation;
public:
    BaseBenchmark(size_t size)
        : InternalItems(size)
        , Permutation(size) 
    {
        for (size_t i = 0; i < size; ++i) {
            InternalItems[i] = i;
            Permutation[i] = i;
        }
        random_shuffle(Permutation.begin(), Permutation.end());
    }

    std::chrono::duration<double> Run(size_t nTests) const {
        using namespace std::chrono;
        vector<int> items = InternalItems;
        high_resolution_clock::time_point startTime = high_resolution_clock::now();
        for (size_t i = 0; i < nTests; ++i) 
            apply_permutation(items, Permutation);
        high_resolution_clock::time_point finTime = high_resolution_clock::now();
        return duration_cast<duration<double>> (finTime - startTime);
    }

    virtual double RunBenchmark(size_t nRuns, size_t testsForRun) const {
        throw domain_error("not implemented");
    }
};


class MinExecTimeBenchmark: public BaseBenchmark {
public:
    MinExecTimeBenchmark(size_t size)
        : BaseBenchmark(size)
    {}

    virtual double RunBenchmark(size_t nRuns, size_t testsForRun) const {
        assert (nRuns > 0);
        vector<double> times(nRuns);
        for (size_t i = 0; i < nRuns; ++i)
            times[i] = Run(testsForRun).count();
        return *min_element(times.begin(), times.end());
    }
};


class MedianTimeBenchmark: public BaseBenchmark {
public:
    MedianTimeBenchmark(size_t size)
        : BaseBenchmark(size)
    {}

    virtual double RunBenchmark(size_t nRuns, size_t testsForRun) const {
        assert (nRuns > 0);
        vector<double> times(nRuns);
        for (size_t i = 0; i < nRuns; ++i)
            times[i] = Run(testsForRun).count();
        sort(times.begin(), times.end());
        return times.size() % 2 ? times[times.size()/2] : (times[times.size()/2] + times[times.size()/2 - 1])/2;
    }
};


class AverageTimeBenchmark: public BaseBenchmark {
public:
    AverageTimeBenchmark(size_t size) 
        : BaseBenchmark(size)
    {}

    virtual double RunBenchmark(size_t nRuns, size_t testsForRun) const {
        assert (nRuns > 0);
        vector<double> times(nRuns);
        for (size_t i = 0; i < nRuns; ++i)
            times[i] = Run(testsForRun).count();
        return accumulate(times.begin(), times.end(), 0.0) / nRuns;
    }
};


class CuttedAverageTimeBenchmark: public BaseBenchmark {
public:
    CuttedAverageTimeBenchmark(size_t size)
        : BaseBenchmark(size)
    {}

    virtual double RunBenchmark(size_t nRuns, size_t testsForRun) const {
        assert (nRuns > 4);
        vector<double> times(nRuns);
        for (size_t i = 0; i < nRuns; ++i)
            times[i] = Run(testsForRun).count();
        sort(times.begin(), times.end());
        return accumulate(times.begin() + 2, times.end() - 2, 0.0) / (times.size() - 4);
    }
};


pair<double, double> ComputeMeanAndVariance(const vector<double>& values) {
    // Wellford method for accurate calculation of mean and variance for vector 'values'
    // see http://www.johndcook.com/standard_deviation.html for explanations
    assert (values.size() > 1);
    double oldM, newM;
    double oldS, newS;

    oldM = newM = values.front();
    oldS = 0.0;

    for (size_t i = 1; i < values.size(); ++i) {
        const double x = values[i];
        newM = oldM + (x - oldM) / i;
        newS = oldS + (x - oldM) * (x - newM);

        oldM = newM;
        oldS = newS;
    }

    return make_pair(newM, newS / (values.size() - 1));
}


void PrintBenchmarkStatistics(const string& benchname, const BaseBenchmark* benchmark, size_t benchRuns, size_t nRuns, size_t testsForRun) {
    vector<double> values;
    for (size_t i = 0; i < benchRuns; ++i) 
        values.push_back(benchmark->RunBenchmark(nRuns, testsForRun));
    pair<double, double> stat = ComputeMeanAndVariance(values);
    cout << "benchmark \"" << benchname << "\"\tmean: " << stat.first << "; variance: " << sqrt(stat.second) 
         << " (" << 100.0 * sqrt(stat.second) / stat.first << "%)" << endl;
}


int main() {
    const size_t ITEMS = 100;
    const size_t WARMUP_TESTS = 10000000;
    const size_t BENCH_TESTS = 100000;

    const size_t BENCH_RUNS = 100;
    const size_t VARIANCE_CALC_RUNS = 10;
    try {
        {
            //run warmup for 1 mln cycles
            BaseBenchmark benchmark(ITEMS);
            cout << "warming up time: " << benchmark.Run(WARMUP_TESTS).count() << endl;
        }

        //calculate average and variance for different benchmark types
        //less variance is better
        PrintBenchmarkStatistics("minexectime", unique_ptr<BaseBenchmark>(new MinExecTimeBenchmark(ITEMS)).get(), VARIANCE_CALC_RUNS, BENCH_RUNS, BENCH_TESTS);
        PrintBenchmarkStatistics("mediantime", unique_ptr<BaseBenchmark>(new MedianTimeBenchmark(ITEMS)).get(), VARIANCE_CALC_RUNS, BENCH_RUNS, BENCH_TESTS);
        PrintBenchmarkStatistics("averagetime", unique_ptr<BaseBenchmark>(new AverageTimeBenchmark(ITEMS)).get(), VARIANCE_CALC_RUNS, BENCH_RUNS, BENCH_TESTS);
        PrintBenchmarkStatistics("cutavgtime", unique_ptr<BaseBenchmark>(new CuttedAverageTimeBenchmark(ITEMS)).get(), VARIANCE_CALC_RUNS, BENCH_RUNS, BENCH_TESTS);
    } catch (const exception& e) {
        cerr << "exception occured: " << e.what() << endl;
    }
    return 0;
}
