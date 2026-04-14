function [V4d, info] = convert_dicom_to_nifti_4d(dcmList, outPath)
    ensure_dir_local(fileparts(outPath));
    V4d = dicom_series_to_4d(dcmList);
    if ndims(V4d) == 3
        V4d = reshape(V4d, size(V4d, 1), size(V4d, 2), size(V4d, 3), 1);
    end
    assert(ndims(V4d) == 4, 'DICOM 转换后功能像应为 4D。');
    niftiwrite(single(V4d), outPath, 'Compressed', false);
    info = niftiinfo(outPath);
    V4d = single(niftiread(info));
    if ndims(V4d) == 3
        V4d = reshape(V4d, size(V4d, 1), size(V4d, 2), size(V4d, 3), 1);
    end
    assert(ndims(V4d) == 4, 'NIfTI 回读后的功能像应为 4D。');
end

function ensure_dir_local(p)
    if ~isfolder(p)
        mkdir(p);
    end
end
