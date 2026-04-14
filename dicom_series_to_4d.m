function V4d = dicom_series_to_4d(dcmList)
    assert(~isempty(dcmList), 'DICOM 列表为空。');
    files = fullfile({dcmList.folder}, {dcmList.name});
    nFiles = numel(files);
    hdr = cell(nFiles, 1);
    inst = zeros(nFiles, 1);
    slicePos = nan(nFiles, 1);
    timeKey = nan(nFiles, 1);

    for i = 1:nFiles
        h = dicominfo(files{i});
        hdr{i} = h;
        inst(i) = read_instance_number(h, i);
        slicePos(i) = read_slice_position(h);
        timeKey(i) = read_time_key(h);
    end

    [~, ord] = sort(inst);
    files = files(ord);
    hdr = hdr(ord);
    slicePos = slicePos(ord);
    timeKey = timeKey(ord);

    validTime = ~isnan(timeKey);
    if any(validTime)
        tVals = unique(timeKey(validTime));
        if numel(tVals) > 1
            [V4d, ok] = build_from_time_groups(files, slicePos, timeKey, tVals);
            if ok
                if ndims(V4d) == 3
                    V4d = reshape(V4d, size(V4d, 1), size(V4d, 2), size(V4d, 3), 1);
                end
                return;
            end
        end
    end

    nZ = infer_n_slices(slicePos, hdr{1});
    assert(mod(nFiles, nZ) == 0, ...
        '无法从 DICOM 序列稳定推断 4D 结构：文件数(%d)不能被每时相切片数(%d)整除。', nFiles, nZ);
    nT = nFiles / nZ;

    s = dicomread(files{1});
    V4d = zeros([size(s), nZ, nT], 'single');
    for t = 1:nT
        idx = (t-1)*nZ + (1:nZ);
        if all(~isnan(slicePos(idx)))
            [~, o] = sort(slicePos(idx));
            idx = idx(o);
        end
        for z = 1:nZ
            V4d(:, :, z, t) = single(dicomread(files{idx(z)}));
        end
    end

    if ndims(V4d) == 3
        V4d = reshape(V4d, size(V4d, 1), size(V4d, 2), size(V4d, 3), 1);
    end
end

function v = read_instance_number(h, fallback)
    if isfield(h, 'InstanceNumber') && ~isempty(h.InstanceNumber)
        v = double(h.InstanceNumber);
    else
        v = fallback;
    end
end

function z = read_slice_position(h)
    z = nan;
    if isfield(h, 'ImagePositionPatient') && numel(h.ImagePositionPatient) >= 3
        z = double(h.ImagePositionPatient(3));
        return;
    end
    if isfield(h, 'SliceLocation') && ~isempty(h.SliceLocation)
        z = double(h.SliceLocation);
    end
end

function t = read_time_key(h)
    t = nan;
    if isfield(h, 'TemporalPositionIdentifier') && ~isempty(h.TemporalPositionIdentifier)
        t = double(h.TemporalPositionIdentifier);
        return;
    end
    if isfield(h, 'AcquisitionNumber') && ~isempty(h.AcquisitionNumber)
        t = double(h.AcquisitionNumber);
        return;
    end
    if isfield(h, 'TriggerTime') && ~isempty(h.TriggerTime)
        t = double(h.TriggerTime);
        return;
    end
    if isfield(h, 'AcquisitionTime') && ~isempty(h.AcquisitionTime)
        t = parse_dicom_hhmmss(h.AcquisitionTime);
        if ~isnan(t), return; end
    end
    if isfield(h, 'ContentTime') && ~isempty(h.ContentTime)
        t = parse_dicom_hhmmss(h.ContentTime);
    end
end

function sec = parse_dicom_hhmmss(v)
    sec = nan;
    if isnumeric(v)
        sec = double(v);
        return;
    end
    s = char(string(v));
    s = regexprep(s, '[^\d\.]', '');
    if isempty(s)
        return;
    end
    p = regexp(s, '^(\d{2})(\d{2})(\d{2}(?:\.\d+)?)$', 'tokens', 'once');
    if isempty(p)
        n = str2double(s);
        if ~isnan(n), sec = n; end
        return;
    end
    hh = str2double(p{1});
    mm = str2double(p{2});
    ss = str2double(p{3});
    if ~isnan(hh) && ~isnan(mm) && ~isnan(ss)
        sec = hh * 3600 + mm * 60 + ss;
    end
end

function [V4d, ok] = build_from_time_groups(files, slicePos, timeKey, tVals)
    ok = false;
    V4d = [];
    nT = numel(tVals);
    idxByT = cell(nT, 1);
    nZ = zeros(nT, 1);
    for t = 1:nT
        idx = find(timeKey == tVals(t));
        if isempty(idx)
            return;
        end
        if all(~isnan(slicePos(idx)))
            [~, o] = sort(slicePos(idx));
            idx = idx(o);
        end
        idxByT{t} = idx;
        nZ(t) = numel(idx);
    end
    if any(nZ ~= nZ(1))
        return;
    end
    s = dicomread(files{idxByT{1}(1)});
    V4d = zeros([size(s), nZ(1), nT], 'single');
    for t = 1:nT
        idx = idxByT{t};
        for z = 1:nZ(1)
            V4d(:, :, z, t) = single(dicomread(files{idx(z)}));
        end
    end
    ok = true;
end

function nZ = infer_n_slices(slicePos, firstHdr)
    nZ = nan;
    valid = ~isnan(slicePos);
    if any(valid)
        sp = slicePos(valid);
        sp = round(sp * 1e3) / 1e3;
        u = unique(sp);
        if numel(u) > 1
            nZ = numel(u);
        end
    end
    if isnan(nZ) || nZ <= 0
        if isfield(firstHdr, 'NumberOfSlices') && ~isempty(firstHdr.NumberOfSlices)
            nZ = double(firstHdr.NumberOfSlices);
        elseif isfield(firstHdr, 'ImagesInAcquisition') && ~isempty(firstHdr.ImagesInAcquisition)
            nZ = double(firstHdr.ImagesInAcquisition);
        else
            nZ = 1;
        end
    end
end
