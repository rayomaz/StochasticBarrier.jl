function sb_record(system_name, cls, tbl, degree, time_s, eta_val, beta_val, Ps, status)
% Append a benchmark result row to $SB_RESULTS_DIR/sostools_results.csv.
%
% Columns:
%   tool,table,class,system,barrier_type,pwc_method,degree,num_partitions,
%   time_s,eta,beta,Ps,status

    out_dir = getenv('SB_RESULTS_DIR');
    if isempty(out_dir)
        out_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    path = fullfile(out_dir, 'sostools_results.csv');
    write_header = ~exist(path, 'file') || dir_filesize(path) == 0;

    fid = fopen(path, 'a');
    if fid == -1
        warning('sb_record: could not open %s for append', path);
        return
    end
    cleaner = onCleanup(@() fclose(fid));

    if write_header
        fprintf(fid, ['tool,table,class,system,barrier_type,pwc_method,' ...
                      'degree,num_partitions,time_s,eta,beta,Ps,status\n']);
    end

    if strcmp(status, 'OK')
        fprintf(fid, 'sostools,%d,%s,%s,SOS,,%d,NA,%.6f,%.6g,%.6g,%.6g,OK\n', ...
                tbl, cls, system_name, degree, time_s, eta_val, beta_val, Ps);
    else
        fprintf(fid, 'sostools,%d,%s,%s,SOS,,%d,NA,%.6f,,,,FAILED\n', ...
                tbl, cls, system_name, degree, time_s);
    end
end

function sz = dir_filesize(path)
    info = dir(path);
    if isempty(info)
        sz = 0;
    else
        sz = info(1).bytes;
    end
end
