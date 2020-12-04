# Basic Usage
1. Copy wherever you want the ```makefile``` and the desired driver (```ed_<modelname>.f90```)
2. Customize to your needs the upper part of the ```makefile```. This amount to setting:
    - The driver name *without* ```.f90``` extension, e.g. ```EXE=ed_bhz_2d```    
    - The fortran compiler [```FC```] of choice; ```ifort```, ```gnu``` (> v 4.7) or ```mpif90``` are supported   
    - The platform of choice; either ```gnu``` or ```intel```   
    - The target directory for the compiled executable [```DIREXE```]; Remember to add it to the ```PATH```*    
3. Open a terminal therein and run ```$ make ```: the executable will be built and placed in ```DIREXE```
4. Open in your working directory a terminal and run ```ed_<modelname>```. 
    - A default input file will be dumped therein as ```used.input<MODELNAME>.conf```
5. Edit to your needs the input file with your text editor of choice and save it as ```input<MODELNAME>.conf```
6. Then you can run again the executable and start crunching numbers!

\* Append ```export PATH=<DIREXE>:$PATH``` to your ```$HOME/.bashrc``` file

# Output Files
- The ```*.restart``` files are useful if you want to run a new calculation initializing the dmft bath according to the final solution of the previous one. Just copypaste all the ```*.restart``` files in the directory of the new calculation, together with the (updated) input file and run: the code will find and use them. 
- The presence of a ```ERROR.README``` file signals that the dmft loop has not converged: number of iterations exceeded its *maximum value* [```NLOOP```], without satisfying the *convergence condition* on the bath, as defined by the ```DMFT_ERROR```; both parameters are specified by the user in the input file. Be aware if it.
