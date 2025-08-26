# Wireless-Communications  
**Project-1---Fading-Channel-Simulator**  
本次作業旨在模擬無線通訊中的 Rayleigh 衰落通道，探討移動環境下信號 因多徑傳播所產生的衰落特性。模擬採用 Jakes 模型，產生符合 U 型 Doppler  頻譜的衰落信號，並透過複數高斯雜訊產生 I/Q 分支，經 FFT 與 Doppler濾波 器處理後轉換回時域以獲得衰落通道。作業進一步比較不同的最大 Doppler頻率 對通道自相關、功率分佈與統計行為的影響。    
**Project 2 - Adaptive Equalizer**  
是透過建立 Rayleigh 衰落加上 AWGN 雜訊的通道模型，設計並比較 LMS 與 RLS 兩種自適應等化器的性能。LMS 演算法以最小均方誤差為目標，結構簡單，收斂速度取決於步長參數； RLS 則藉由遞迴最小平方法快速收斂，但計算複雜度較高。透過 MSE 曲線、steady-state misadjustment、以及 BER vs. SNR 分析來進行兩種演算法的性能評估。     
**Project 3 Diversity vs. Multiplexing in a 2×2 MIMO Link**
比較兩種典型的 MIMO 傳輸模式：Alamouti 空間時域分集（STBC) 與空間多工（Spatial Multiplexing） 並在在 flat Rayleigh fading 通道下分析其 BER 表現與頻譜效率差異。 
