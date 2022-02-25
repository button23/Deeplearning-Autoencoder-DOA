

[vec]=DL_SCM_2_vec(h_Rxx,M0); % noisy SCM to vector
dataName = "vec";
rootPath = 'C:\Users\HYPC300\OneDrive - 한양대학교\GitHub\Deeplearning-Autoencoder-DOA\data\vec';
save(rootPath, dataName) % Save noisy data.
