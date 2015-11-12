disp('Executing time vectorized code: ')
tic; rec=fdtd2ds(0); toc;
disp('------')
disp('Executing time code modified with 1/dy => dy to avoid divisions: ')
tic; rec=fdtd2ds_oneOver(0); toc;
disp('------')
disp('Executing time loop code')
tic; rec=fdtd2ds_loop(0); toc;
disp('------')
exit;
