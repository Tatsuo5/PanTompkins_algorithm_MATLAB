# Pan-Tompkins Algorithm Implementation in MATLAB

A MATLAB implementation of the Pan-Tompkins algorithm for QRS detection in ECG signals.

## Features

- Detects QRS complexes in ECG signals using the Pan-Tompkins algorithm.
- Includes preprocessing steps such as bandpass filtering, differentiation, squaring, and moving average smoothing.
- Executes re-search and T-wave identification as needed, and performs adaptive threshold updates after each peak detection.
- Excludes abnormal peaks from threshold updates based on user-defined criteria.
- Accurately identifies QRS complexes by distinguishing R-waves from T-waves.
- Provides customizable visualization of R-peaks, adaptive thresholds, and T-wave exclusion.

## Usage

### Input

- **Sampling frequency (fs)**: The frequency at which the ECG signal is sampled, specified in Hz.
- **Raw data (data)**: A matrix containing the ECG signal and possibly other data.
- **ECG column index (ecg_column)**: The column number in the raw data matrix where the ECG signal is located.
- **Threshold multiplier (k)**: A user-defined value to determine how abnormal peaks are excluded from threshold updates.  
  Recommended values are approximately 50–100. If you do not wish to exclude abnormal peaks from threshold updates, set `k` to a very high value (e.g., around 1 million).  
  This parameter specifies that when a peak's amplitude exceeds the current threshold by a factor of `k`, it will still be detected as a peak, but its amplitude will not be used to update the adaptive threshold.
- **Plot flag (plot_flag)**: A flag to enable or disable plot output. Set to 1 for visualization or 0 to disable plotting.

### Output

The `Pan_Tompkins` function returns the following outputs:

- **peaks_idx**: Detected R-wave indices.
- **peaks_data**: Amplitudes (magnitudes) of the detected R-waves.
- **smoothed_data**: Preprocessed ECG signal.
- **THRESHOLD_I1**: Moving threshold signal (first threshold, with the second threshold being half of this).


## Algorithm Details

The Pan-Tompkins algorithm processes ECG signals through the following steps:

1. **Preprocessing**:
   - **Bandpass filtering**: Removes noise and baseline wander by retaining frequencies in the range of 5-15 Hz.
   - **Differentiation filter**: Highlights the slope of QRS complexes using the differentiation filter proposed by Pan and Tompkins:

     $$
     H(z) = \frac{1}{8} \left( -z^{-2} - 2z^{-1} + 2z^{1} + z^{2} \right)
     $$

   - **Squaring**: Amplifies high-frequency components and suppresses negative amplitudes to emphasize QRS peaks.
   - **Moving average smoothing**: Extracts the signal envelope to further highlight QRS peaks and reduce noise.

2. **Peak Detection**:
   - Identifies R-peaks using adaptive thresholds.

   **Note:** The initial detected R-peak data may contain false detections due to the adaptation of the threshold. It is recommended to exclude the first detected R-peak(s) when analyzing the results.

3. **Re-search**:
   - If a QRS complex is not detected within 166% of the average RR interval (calculated from the most recent 8 detected R-peaks), the algorithm searches for the largest peak in that interval.
   - The re-search process uses a secondary threshold (THRESHOLD_I2) to increase sensitivity while ensuring false detections are minimized.
   - This step ensures robustness in detecting QRS complexes, even in cases of irregular heartbeats or missed peaks.

4. **T-Wave Identification**:
   - When an RR interval is less than 360 ms and greater than the 200 ms latency period, the algorithm checks whether the current peak is a T-wave or a QRS complex:
     - The maximum slope of the peak is compared to the maximum slope of the preceding QRS complex.
     - If the maximum slope is less than half of the preceding QRS maximum slope, it is classified as a T-wave.
     - Otherwise, it is classified as a QRS complex.

5. **Threshold Updates**:
   - **Initial Threshold Calculation**:
     - During the learning phase (first 2 seconds of the signal), the algorithm calculates the initial thresholds:
       - **SPKI (signal peak indicator)**: The mean of peaks detected in the learning phase.
       - **NPKI (noise peak indicator)**: The mean of peaks below the 50th percentile in the learning phase.
   - **Threshold Update Equations**:
     - After every detected R-peak, the thresholds are updated as follows:

       THRESHOLD\_I1 = NPKI + 0.25 (SPKI - NPKI)

       THRESHOLD\_I2 = 0.5 THRESHOLD\_I1

     - Where:
       - **SPKI** is updated with the amplitude of detected signal peaks.
       - **NPKI** is updated with the amplitude of detected noise peaks.

6. **Enhancements**:
   - Excludes abnormal peaks from threshold updates based on user-defined criteria.
   - Provides customizable visualization of results, including detected R-peaks and adaptive thresholds.

## Sample Code

The sample code demonstrating the usage of the `Pan_Tompkins` function and block segmentation of ECG data is available in the repository under the file `sample_code.m`.

Additionally, a sample ECG dataset (`sample_ecg.csv`) is provided in the repository. This dataset can be used to test the algorithm and visualize its functionality.

### Highlights:
- How to run the `Pan_Tompkins` function on sample ECG data.
- Block segmentation of intervals and RR interval calculation.
- Examples of setting parameters like sampling frequency, signal column, and thresholds.

Refer to the provided `sample_code.m` file and `sample_ecg.csv` dataset for detailed implementation and testing.

## License

This project is licensed under the MIT License. You are free to use, modify, and distribute this code, provided proper attribution is given. See the LICENSE file in the repository for more details.

## Author

This project was created by **Tatsuo Hirata**.  
For any questions or feedback, please use the Issues section of this repository.

## Disclaimer

This software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, or non-infringement. The authors are not responsible for any damage or loss caused by the use of this software. Users are encouraged to validate results independently.
