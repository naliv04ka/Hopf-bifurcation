#ifndef FFT_ANALYZER_H
#define FFT_ANALYZER_H

#include "common.h"

class FFTAnalyzer {
public:
    FFTAnalyzer();
    
    // Add data point
    void addDataPoint(Real time, Real value);
    
    // Compute FFT
    void computeFFT();
    
    // Get dominant frequency and amplitude
    Real getDominantFrequency() const { return dominantFreq_; }
    Real getDominantAmplitude() const { return dominantAmp_; }
    
    // Get frequency spectrum
    const std::vector<Real>& getFrequencies() const { return frequencies_; }
    const std::vector<Real>& getMagnitudes() const { return magnitudes_; }
    
    // Save data
    void saveTimeSeries(const std::string& filename) const;
    void saveSpectrum(const std::string& filename) const;
    
    // Clear data
    void clear();
    
private:
    std::vector<Real> times_;
    std::vector<Real> values_;
    std::vector<Real> frequencies_;
    std::vector<Real> magnitudes_;
    
    Real dominantFreq_;
    Real dominantAmp_;
    
    // FFT implementation
    void fft(std::vector<std::complex<Real>>& data);
    void fftRecursive(std::vector<std::complex<Real>>& data);
    
    // Find nearest power of 2
    Index nearestPowerOf2(Index n);
};

#endif // FFT_ANALYZER_H
