function data_pca = performPCA(data)
    [~, score, ~] = pca(data);
    data_pca = score(:, 1);
end