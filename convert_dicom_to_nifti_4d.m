function [V4d, info] = convert_dicom_to_nifti_4d(dcmList, outPath)
    ensure_dir_local(fileparts(outPath));
    V4d = ensure_4d_local(dicom_series_to_4d(dcmList));
    assert(ndims(V4d) == 4, 'DICOM 转换后功能像应为 4D。');
    niftiwrite(single(V4d), outPath, 'Compressed', false);
    info = niftiinfo(outPath);
    V4d = ensure_4d_local(single(niftiread(info)));
    assert(ndims(V4d) == 4, 'NIfTI 回读后的功能像应为 4D。');
end

function ensure_dir_local(p)
    if ~isfolder(p)
        mkdir(p);
    end
end

function V = ensure_4d_local(V)
    if ndims(V) == 3
        V = reshape(V, size(V, 1), size(V, 2), size(V, 3), 1);
    end
end
