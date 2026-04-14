function out = read_ctex_pdf_text(pdfPath, outTxtPath)
% READ_CTEX_PDF_TEXT
% -------------------------------------------------------------------------
% 针对中文 ctex PDF 的稳健文本提取工具。
% 提取顺序：
% 1) MATLAB extractFileText
% 2) 系统 pdftotext -enc UTF-8 -layout
% 3) Python + pypdf
% -------------------------------------------------------------------------

    if nargin < 1 || isempty(pdfPath)
        repoRoot = fileparts(mfilename('fullpath'));
        pdfPath = fullfile(repoRoot, 'MRILAB3.pdf');
    end
    if nargin < 2
        outTxtPath = '';
    end

    assert(exist(pdfPath, 'file') == 2, 'PDF not found: %s', pdfPath);

    out = struct('success', false, 'method', '', 'text', '');

    % 方法1：MATLAB 内置解析
    try
        t = extractFileText(pdfPath);
        if is_usable_text(t)
            out.success = true;
            out.method = 'extractFileText';
            out.text = t;
        end
    catch
    end

    % 方法2：系统 pdftotext
    if ~out.success
        [ok, t] = extract_with_pdftotext(pdfPath);
        if ok && is_usable_text(t)
            out.success = true;
            out.method = 'pdftotext';
            out.text = t;
        end
    end

    % 方法3：Python pypdf
    if ~out.success
        [ok, t] = extract_with_pypdf(pdfPath);
        if ok && is_usable_text(t)
            out.success = true;
            out.method = 'python_pypdf';
            out.text = t;
        end
    end

    if ~out.success
        error('PDF text extraction failed for: %s', pdfPath);
    end

    if ~isempty(outTxtPath)
        fid = fopen(outTxtPath, 'w');
        assert(fid ~= -1, 'Cannot write file: %s', outTxtPath);
        fwrite(fid, out.text, 'char');
        fclose(fid);
    end
end

function tf = is_usable_text(t)
    tf = ischar(t) || isstring(t);
    if ~tf, return; end
    t = char(t);
    tf = strlength(string(strtrim(t))) > 100;
end

function [ok, textOut] = extract_with_pdftotext(pdfPath)
    cmd = sprintf('pdftotext -enc UTF-8 -layout "%s" -', pdfPath);
    [status, out] = system(cmd);
    ok = (status == 0);
    textOut = out;
end

function [ok, textOut] = extract_with_pypdf(pdfPath)
    py = [ ...
        "import pypdf,sys;" ...
        + "r=pypdf.PdfReader(r'" + string(pdfPath) + "');" ...
        + "sys.stdout.write('\\n'.join([(p.extract_text() or '') for p in r.pages]))" ...
    ];
    cmd = sprintf('python -c "%s"', escape_double_quotes(py));
    [status, out] = system(cmd);
    ok = (status == 0);
    textOut = out;
end

function s = escape_double_quotes(s)
    s = strrep(char(s), '"', '\"');
end
