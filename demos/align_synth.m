run('../set_path.m');

[img_noise, v_ref, ref, img] = generate_synth_ERP(35);
% imagesc(img);

alpha = [5, 4, 3, 1];
sigma = [10, 8, 5, 2];
sigma_l = [15, 15, 10, 8];

% alpha = 1;
% sigma = 2;
% sigma_l = 8;

img_reg = img_noise;
ref_init = mean(img_reg, 1);
hold on
plot(ref);
plot(ref_init);

img_reg_vid = zeros([size(img_reg), length(alpha) + 1]);
img_reg_vid(:, :, 1) = img_reg;
v = zeros(size(img_reg));
for i = 1:length(alpha)
    reg_filt = imgaussfiltaniso(img_reg, sigma(i));

    if i == 1
        ref_tmp = mat2gray(mean(reg_filt(1:50, :), 1));
    else
        ref_tmp = mat2gray(mean(reg_filt, 1));
    end

    [~, v_tmp] = align_lines(...
        reg_filt, ...
        ref_tmp, ...
        'sigma', sigma_l(i), ...
        'iterations', 50, ...
        'alpha', alpha(i), ...
        'a_smooth', 1);

    img_reg = horiz_alignment(img_noise, v_tmp + v);
    img_reg_vid(:, :, i + 1) = img_reg;
    v = v + v_tmp;
end
plot(mean(img_reg, 1))

figure
subplot(1, 2, 1);
imagesc(v);
subplot(1, 2, 2);
imagesc(v_ref);