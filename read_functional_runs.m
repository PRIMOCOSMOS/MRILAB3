function [runs, infos] = read_functional_runs(funcDir, outDir)
    ensure_dir_local(outDir);
    runs = {};
    infos = {};

    niiTop = dir(fullfile(funcDir, '*.nii*'));
    if ~isempty(niiTop)
        p = fullfile(niiTop(1).folder, niiTop(1).name);
        info = niftiinfo(p);
        V = single(niftiread(info));
        assert(ndims(V) == 4, '功能像 NIfTI 必须为 4D。');
        runs{1} = V; infos{1} = info; %#ok<AGROW>
        write_nifti_4d_local(V, info, fullfile(outDir, 'run01_raw.nii'));
        return;
    end

    d = dir(funcDir);
    d = d([d.isdir]);
    d = d(~startsWith({d.name}, '.'));
    if isempty(d)
        dcmTop = dir(fullfile(funcDir, '**', '*.dcm'));
        if ~isempty(dcmTop)
            outPath = fullfile(outDir, 'run01_raw.nii');
            [V, info] = convert_dicom_to_nifti_4d(dcmTop, outPath);
            runs{1} = V; infos{1} = info; %#ok<AGROW>
            return;
        end
    end
    for i = 1:numel(d)
        runDir = fullfile(funcDir, d(i).name);
        nii = dir(fullfile(runDir, '*.nii*'));
        if ~isempty(nii)
            p = fullfile(nii(1).folder, nii(1).name);
            info = niftiinfo(p);
            V = single(niftiread(info));
            assert(ndims(V) == 4, 'run %s 需为4D。', d(i).name);
            write_nifti_4d_local(V, info, fullfile(outDir, sprintf('run%02d_raw.nii', i)));
        else
            dcm = dir(fullfile(runDir, '**', '*.dcm'));
            assert(~isempty(dcm), 'run %s 目录中无 NIfTI 或 DICOM。', d(i).name);
            outPath = fullfile(outDir, sprintf('run%02d_raw.nii', i));
            [V, info] = convert_dicom_to_nifti_4d(dcm, outPath);
        end
        runs{end + 1} = V; %#ok<AGROW>
        infos{end + 1} = info; %#ok<AGROW>
    end
end

function ensure_dir_local(p)
    if ~isfolder(p)
        mkdir(p);
    end
end

function write_nifti_4d_local(vol4d, refInfo, pathOut)
    ensure_dir_local(fileparts(pathOut));
    info = sanitize_nifti_info_for_write_local(refInfo);
    info.ImageSize = size(vol4d);
    info.Datatype = 'single';
    safe_niftiwrite_local(single(vol4d), pathOut, info);
end

function infoOut = sanitize_nifti_info_for_write_local(infoIn)
    infoOut = struct();
    if ~isstruct(infoIn)
        return;
    end
    keep = {'ImageSize','PixelDimensions','Datatype','SpaceUnits','TimeUnits', ...
            'AdditiveOffset','MultiplicativeScaling','Transform','Qfactor'};
    for i = 1:numel(keep)
        k = keep{i};
        if isfield(infoIn, k)
            infoOut.(k) = infoIn.(k);
        end
    end
end

function safe_niftiwrite_local(vol, pathOut, info)
    infoTry = info;
    minStripAttempts = 4;
    extraStripAttempts = 2;
    maxStripAttempts = max(minStripAttempts, numel(fieldnames(infoTry)) + extraStripAttempts);
    for k = 1:maxStripAttempts
        try
            niftiwrite(vol, pathOut, infoTry, 'Compressed', false);
            return;
        catch ME
            badField = extract_unknown_info_field_local(ME.message);
            if ~isempty(badField) && isfield(infoTry, badField)
                infoTry = rmfield(infoTry, badField);
                continue;
            end
            break;
        end
    end
    warning('niftiwrite:InfoFallback', ...
        'Failed to keep compatible NIfTI header for %s; fallback to write without Info.', pathOut);
    niftiwrite(vol, pathOut, 'Compressed', false);
end

function fieldName = extract_unknown_info_field_local(msg)
    fieldName = '';
    quoteClass = '"“”''';
    cnUnknownFieldPattern = ['无法识别的字段名称[:：\s]*[', quoteClass, ']?([^', quoteClass, '\s]+)'];
    t = regexp(msg, cnUnknownFieldPattern, 'tokens', 'once');
    if ~isempty(t)
        fieldName = t{1};
        return;
    end
    enUnknownFieldPattern = 'Unrecognized field name\s*''([^'']+)''';
    t = regexp(msg, enUnknownFieldPattern, 'tokens', 'once');
    if ~isempty(t)
        fieldName = t{1};
    end
end
