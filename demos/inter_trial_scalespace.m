run('../set_path.m');

[img_noise, v_ref, ref, img] = generate_synth_ERP(35);


img_reg = img_noise;
ref_init = mean(img_reg, 1);
hold on
plot(ref);
plot(ref_init);

sigma_aniso = 5;
sigma = 10;
alpha = 3;

trial_eta = 0.75;
trial_levels = 10;
[n, m] = size(img_reg);

max_level = warpingDepthY(trial_eta, trial_levels, n);
v = zeros(size(img_reg));


for i = max_level:-1:1
    f1_level = imresize(img_reg, [n * trial_eta^i, m], ...
        'Antialiasing', true);

    reg_filt = imgaussfiltaniso(f1_level, sigma);
    ref_tmp = mean(img_reg, 1);

    [~, v_tmp] = align_lines(...
        reg_filt, ...
        ref_tmp, ...
        'sigma', 5, ...
        'iterations', 50, ...
        'alpha', alpha, ...
        'a_smooth', 1);

    v_tmp = imresize(v_tmp, [n, m]);
    img_reg = horiz_alignment(img_noise, v_tmp + v);

    v = v_tmp + v;
end

subplot(1, 3, 2);
imagesc(img_reg);
subplot(1, 3, 1);
imagesc(img_noise);
subplot(1, 3, 3);
plot(squeeze(mean(img_noise, 1)));
hold on
plot(squeeze(mean(img_reg, 1)));
hold off


function d = warpingDepthY(eta, levels, m)
    warpingdepth = 0;
    d = warpingdepth;

    for i = 1:levels
        warpingdepth = warpingdepth + 1;
        m = m * eta;
        if (round(m) < 10	)
            break;
        end
        d = warpingdepth;
    end
end