%% === PANIC-SAVE of `out` ===
assert(exist('out','var')==1, '`out` not found in workspace.');

ts = datestr(now,'yyyymmdd_HHMMSS');
tag = sprintf('mask_sweep_backup_%s', ts);
outdir = pwd;  % change if you want a subfolder
fprintf('Saving panic backups with tag: %s\n', tag);

% 1) Full MAT backup (best fidelity; keeps LTI objects)
matfile = fullfile(outdir, [tag '.mat']);
try
  save(matfile, 'out', '-v7.3');
  fprintf('Saved MAT: %s\n', matfile);
catch ME
  warning('MAT save failed: %s', ME.message);
end

% 2) Masks table CSV (one row per mask)
try
  Nm = size(out.masks_r,1);
  q  = size(out.masks_r,2);
  p  = size(out.masks_c,2);

  % to human-readable strings like "10101"
  row2str = @(row) string(char('0' + row(:)'));     % logical/0-1 -> "0101..."
  Rstr = strings(Nm,1); Cstr = strings(Nm,1);
  for k=1:Nm
    Rstr(k) = row2str(out.masks_r(k,:));
    Cstr(k) = row2str(out.masks_c(k,:));
  end

  T_masks = table( (1:Nm).', Rstr, Cstr, ...
                   out.wstar(:), out.sigmax(:), ...
                   logical(out.self_destab(:)), ...
                   out.self_min_sigma(:), ...
                   'VariableNames', {'mask_id','rmask','cmask','wstar','sigmax','self_destab','self_min_sigma'});

  masks_csv = fullfile(outdir, [tag '_masks.csv']);
  writetable(T_masks, masks_csv);
  fprintf('Saved masks CSV: %s  (Nm=%d, q=%d, p=%d)\n', masks_csv, Nm, q, p);
catch ME
  warning('Masks CSV save failed: %s', ME.message);
end

% 3) Cross-check CSV (all pairs, stable+unstable)
try
  destab = logical(out.destab);
  minsig = out.min_sigma;
  [J,I] = ndgrid(1:size(destab,1), 1:size(destab,2));  % J: M_j, I: Delta_i
  T_edges = table( I(:), J(:), destab(:), minsig(:), ...
                   'VariableNames', {'delta_id','M_id','destabilizes','min_sigma'} );

  edges_csv = fullfile(outdir, [tag '_destab.csv']);
  writetable(T_edges, edges_csv);
  fprintf('Saved edges CSV: %s  (pairs=%d)\n', edges_csv, numel(destab));
catch ME
  warning('Edges CSV save failed: %s', ME.message);
end

% 4) (Optional) dense matrices as CSVs: destab (0/1) and min_sigma
try
  destab_csv = fullfile(outdir, [tag '_destab_matrix.csv']);
  writematrix(double(out.destab), destab_csv);
  fprintf('Saved destab matrix CSV: %s\n', destab_csv);
catch ME
  warning('destab matrix CSV save failed: %s', ME.message);
end

try
  minsig_csv = fullfile(outdir, [tag '_min_sigma_matrix.csv']);
  writematrix(out.min_sigma, minsig_csv);
  fprintf('Saved min_sigma matrix CSV: %s\n', minsig_csv);
catch ME
  warning('min_sigma matrix CSV save failed: %s', ME.message);
end

% 5) (Optional) compact JSON (metadata only; not LTI objects)
try
  meta = struct;
  meta.timestamp = ts;
  meta.Nm = size(out.masks_r,1);
  meta.q  = size(out.masks_r,2);
  meta.p  = size(out.masks_c,2);
  meta.masks_csv = exist('masks_csv','var')*string(masks_csv) + (~exist('masks_csv','var'))*missing;
  meta.edges_csv = exist('edges_csv','var')*string(edges_csv) + (~exist('edges_csv','var'))*missing;
  meta.destab_fraction = mean(out.destab(:));
  jsonfile = fullfile(outdir, [tag '_meta.json']);
  fid = fopen(jsonfile, 'w'); fwrite(fid, jsonencode(meta, 'PrettyPrint', true)); fclose(fid);
  fprintf('Saved meta JSON: %s\n', jsonfile);
catch ME
  warning('JSON save failed: %s', ME.message);
end

fprintf('Panic-save complete. Keep the MAT file for full fidelity.\n');
