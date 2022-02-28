[vec]=DL_SCM_2_vec(h_Rxx,M0); % noisy SCM to vector
dataName = "vec";
basePath = '8_antenna_500_samples_100_snapshots';
% rootPath = 'C:\Users\HYPC300\OneDrive - 한양대학교\GitHub\Deeplearning-Autoencoder-DOA\data\vec';
rootPath = fullfile('/home/hymc/[0]Github/Deeplearning-Autoencoder-DOA/[2]Tensorflow/[1]DAE_DOA_implementation/vec');
save(rootPath, dataName) % Save noisy data.
