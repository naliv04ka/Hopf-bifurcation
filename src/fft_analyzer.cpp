#include "fft_analyzer.h"
#include <algorithm>
#include <fstream>

FFTAnalyzer::FFTAnalyzer()
    : dominantFreq_(0.0), dominantAmp_(0.0)
{
}

void FFTAnalyzer::addDataPoint(Real time, Real value) {
    times_.push_back(time);
    values_.push_back(value);
}

void FFTAnalyzer::clear() {
    times_.clear();
    values_.clear();
    frequencies_.clear();
    magnitudes_.clear();
    dominantFreq_ = 0.0;
    dominantAmp_ = 0.0;
}

Index FFTAnalyzer::nearestPowerOf2(Index n) {
    Index power = 1;
    while (power < n) {
        power *= 2;
    }
    return power;
}

void FFTAnalyzer::computeFFT() {
    if (values_.empty()) return;
    
    // Prepare data
    Index n = nearestPowerOf2(values_.size());
    std::vector<std::complex<Real>> data(n);
    
    // Remove mean (DC component)
    Real mean = 0.0;
    for (Real val : values_) {
        mean += val;
    }
    mean /= values_.size();
    
    // Copy data and pad with zeros if necessary
    for (Index i = 0; i < n; ++i) {
        if (i < static_cast<Index>(values_.size())) {
            data[i] = std::complex<Real>(values_[i] - mean, 0.0);
        } else {
            data[i] = std::complex<Real>(0.0, 0.0);
        }
    }
    
    // Perform FFT
    fft(data);
    
    // Compute sampling frequency
    Real dt = 0.0;
    if (times_.size() > 1) {
        dt = times_[1] - times_[0]; // Assume uniform sampling
    }
    Real samplingFreq = 1.0 / dt;
    
    // Compute magnitudes and frequencies
    frequencies_.clear();
    magnitudes_.clear();
    
    // Only use first half (Nyquist)
    for (Index i = 0; i < n/2; ++i) {
        Real freq = i * samplingFreq / n;
        Real magnitude = std::abs(data[i]) / n;
        
        frequencies_.push_back(freq);
        magnitudes_.push_back(magnitude);
    }
    
    // Find dominant frequency (skip DC component at i=0)
    dominantAmp_ = 0.0;
    dominantFreq_ = 0.0;
    
    for (Index i = 1; i < static_cast<Index>(magnitudes_.size()); ++i) {
        if (magnitudes_[i] > dominantAmp_) {
            dominantAmp_ = magnitudes_[i];
            dominantFreq_ = frequencies_[i];
        }
    }
}

void FFTAnalyzer::fft(std::vector<std::complex<Real>>& data) {
    Index n = data.size();
    
    if (n <= 1) return;
    
    // Bit-reversal permutation
    Index j = 0;
    for (Index i = 0; i < n-1; ++i) {
        if (i < j) {
            std::swap(data[i], data[j]);
        }
        
        Index m = n / 2;
        while (j >= m) {
            j -= m;
            m /= 2;
        }
        j += m;
    }
    
    // Cooley-Tukey FFT
    for (Index s = 1; s <= std::log2(n); ++s) {
        Index m = 1 << s; // 2^s
        std::complex<Real> wm = std::exp(std::complex<Real>(0, -2.0 * PI / m));
        
        for (Index k = 0; k < n; k += m) {
            std::complex<Real> w(1.0, 0.0);
            
            for (Index j = 0; j < m/2; ++j) {
                std::complex<Real> t = w * data[k + j + m/2];
                std::complex<Real> u = data[k + j];
                
                data[k + j] = u + t;
                data[k + j + m/2] = u - t;
                
                w *= wm;
            }
        }
    }
}

void FFTAnalyzer::saveTimeSeries(const std::string& filename) const {
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    
    file << "# Time\tValue\n";
    for (size_t i = 0; i < times_.size(); ++i) {
        file << times_[i] << "\t" << values_[i] << "\n";
    }
    
    file.close();
}

void FFTAnalyzer::saveSpectrum(const std::string& filename) const {
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    
    file << "# Frequency\tMagnitude\n";
    for (size_t i = 0; i < frequencies_.size(); ++i) {
        file << frequencies_[i] << "\t" << magnitudes_[i] << "\n";
    }
    
    file.close();
}
