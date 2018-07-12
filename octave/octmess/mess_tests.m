%% add source directory to path
addpath(genpath("@CMAKE_SOURCE_DIR@/octave/octmess/src"))


%% call mess_path for autoload
mess_path


%% call tests
test("MESS_matrix")


