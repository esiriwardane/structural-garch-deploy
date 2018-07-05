function lintFiles
  matlabFiles = findMatlabFiles(fullfile('.'));

  checkcode(matlabFiles{:});
end

function matlabFiles = findMatlabFiles(directory)
  matlabFiles = {};

  dirinfo = dir(directory);
  dirinfo(ismember({dirinfo.name}, {'.', '..'})) = [];

  for index = 1:length(dirinfo)
    if dirinfo(index).isdir
      matlabFiles = [ matlabFiles findMatlabFiles(fullfile(directory, dirinfo(index).name)) ]; %#ok
    elseif strcmp(dirinfo(index).name(end - 1:end), '.m')
      matlabFiles = [ matlabFiles fullfile(directory, dirinfo(index).name) ]; %#ok
    end
  end
end