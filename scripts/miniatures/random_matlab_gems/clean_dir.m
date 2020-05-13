function cdir=clean_dir(dir)
bdir=arrayfun(@(x) x.name(1)=='.',dir);
cdir=dir(~bdir);
end