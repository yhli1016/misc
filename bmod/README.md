# bmod

## 1. 简介

在UNIX/Linux系统下做计算，不可避免地会接触到各种环境变量。编译程序时，需要通过环境变量指定编译器、编译选项和依赖库的位置。运行程序时，需要通过环境变量来查找可执行文件和动态链接库。若设置不当，便会导致各种错误（详见附录）。因此管理环境变量便成了一件既重要又棘手的事。

最简单的管理环境变量的方法，就是直接在`$HOME/.bashrc`中添加相应设置。这种方法的缺点也是很明显的：环境变量名易写错，添加或删除设置后要重新登录才能起效，bashrc冗长等等。大型超算中心通常用Environment Modules系统解决这个问题。但目前主流的Environment Modules系统安装和使用都比较复杂，学习成本较高。而个人电脑、自己组装的工作站和小型集群软件环境简单，再用Environment Modules系统牛刀杀鸡之嫌。因此，在这里提供一个环境变量管理程序bmod。程序非常简单，全部由bash语言写成，主体仅有150行左右，但实现了一个Environment Modules系统最基本的功能。下面介绍安装和使用方法。

## 2. 安装

程序下载地址为<https://github.com/yhli1016/misc/tree/master/bmod>。解压后将bmod目录移动到安装路径下面（本文中为`$HOME/soft`）。打开`$HOME/.bashrc`，添加如下设置：

```bash
source $HOME/soft/bmod/init.sh
```

然后执行命令`source $HOME/.bashrc`即可完成安装。

## 3. 使用

### 3.1 通过命令行添加和移除设置

#### 3.1.1 使用预定义的环境变量

与程序编译和运行相关的环境变量众多，不仅难记，还容易拼错。因此，bmod中预定义了常见的环境变量。以安装在`$HOME/soft/lib/fftw-3.3.8`下面的fftw库为例，为使该库能正常工作，我们需要手动输入如下命令：

```bash
fftw_root=$HOME/soft/lib/fftw-3.3.8
export PATH=$fftw_root/bin:$PATH
export LIBRARY_PATH=$fftw_root/lib64:$LIBRARY_PATH
export LD_RUN_PATH=$fftw_root/lib64:$LD_RUN_PATH
export LD_LIBRARY_PATH=$fftw_root/lib64:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$fftw_root/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$fftw_root/include:$CPLUS_INCLUDE_PATH
unset fftw_root
```

而在bmod中，上述设置只需一行命令即可：`set_mod add pkg $HOME/soft/lib/fftw-3.3.8`，移除设置只需`set_mod rm pkg $HOME/soft/lib/fftw-3.3.8`。`set_mod`命令后面第一个参数为`add`或`rm`，指定添加还是删除设置；第二个参数为预定义的环境变量类型；第三个参数为添加到环境变量中，或从环境变量中删除的设置。当预定义类型为`pkg`时，第三个参数无需具体到`bin`、`lib64`或`include`等子目录，bmod会自动搜索和添加。目前bmod中预定义的类型和对应的环境变量为：

- bin
  - PATH
- lib:
  - LIBRARY\_PATH
  - LD\_RUN\_PATH
  - LD\_LIBRARY\_PATH
- inc
  - C\_INCLUDE\_PATH
  - CPLUS\_INCLUDE\_PATH
- py
  - PYTHONPATH
- pkg
  - PATH
  - LIBRARY\_PATH
  - LD\_RUN\_PATH
  - LD\_LIBRARY\_PATH
  - C\_INCLUDE\_PATH
  - CPLUS\_INCLUDE\_PATH

当预定义类型为`bin`、`lib`、`inc`、`py`时，bmod不会自动搜索子目录，所以第三个参数必须具体到`bin`、`lib64`或`include`。上面fftw的例子若用`bin`、`lib`、`inc`等预定义类型改写，对应的命令为：

```bash
set_mod add bin $HOME/soft/lib/fftw-3.3.8/bin
set_mod add lib $HOME/soft/lib/fftw-3.3.8/lib64
set_mod add inc $HOME/soft/lib/fftw-3.3.8/include
```

再看一个xcrysden的例子。程序安装在`$HOME/soft/dft/xcrysden-1.6.2-bin-shared`，可执行文件为`xcrysden`，还有一个动态链接库`libTogl.so.2`。为使这个程序正常运行，需将`$HOME/soft/dft/xcrysden-1.6.2-bin-shared`添加到环境变量`PATH`和`LD_LIBRARY_PATH`，对应操作为：

```bash
set_mod add bin $HOME/soft/dft/xcrysden-1.6.2-bin-shared
set_mod add lib $HOME/soft/dft/xcrysden-1.6.2-bin-shared
```

或者简写为`set_mod add bin+lib $HOME/soft/dft/xcrysden-1.6.2-bin-shared`。

#### 3.1.2 使用自定义的环境变量

预定义的环境变量可以满足大多数情形。若待修改的环境变量没有预定义，可以用`set_env`和`reset_env`命令修改。

假设我们要把`$HOME/soft/lib/abc/def`添加到环境变量`TEST`中，对应操作为`set_env add TEST $HOME/soft/lib/abc/def`，删除时为`set_env rm TEST $HOME/soft/lib/abc/def`。如果我们希望把环境变量重设为一个新的值，而不是把新的值追加到变量中，可以用`reset_env`命令。例如，`reset_env add OMP_NUM_THREADS 4`将OpenMP线程数设置为4，`reset_env rm OMP_NUM_THREADS`将其恢复默认值。

### 3.2 通过脚本添加和移除设置

一般来说，程序在编译和运行过程中会涉及很多库，需要我们输入很多次`set_mod`和`set_env`命令。以siesta为例，加载设置时我们需要输入：

```bash
siesta=$HOME/soft/dft/siesta-v4.1-b4
set_mod add pkg $siesta/Docs/build/flook/0.8.1
set_mod add pkg $siesta/Docs/build/hdf5/1.8.21
set_mod add pkg $siesta/Docs/build/netcdf/4.7.4
set_mod add pkg $siesta/Docs/build/zlib/1.2.11
set_mod add bin $siesta/Obj
unset siesta
reset_env add OMP_NUM_THREADS 1
```

移除设置时需要输入：

```bash
siesta=$HOME/soft/dft/siesta-v4.1-b4
set_mod rm pkg $siesta/Docs/build/flook/0.8.1
set_mod rm pkg $siesta/Docs/build/hdf5/1.8.21
set_mod rm pkg $siesta/Docs/build/netcdf/4.7.4
set_mod rm pkg $siesta/Docs/build/zlib/1.2.11
set_mod rm bin $siesta/Obj
unset siesta
reset_env rm OMP_NUM_THREADS 1
```

而这些命令的区别，仅仅是把`add`换成了`rm`。这样太麻烦了。在bmod中我们可以通过脚本，用同一套命令完成加载和卸载两种操作。

我们在`$HOME/soft/bmod/modules`下面新建一个bash脚本siesta-v4.1-b4.sh，把上面的命令存进去，并把`add`或`rm`换成`$1`：

```bash
siesta=$HOME/soft/dft/siesta-v4.1-b4
set_mod $1 pkg $siesta/Docs/build/flook/0.8.1
set_mod $1 pkg $siesta/Docs/build/hdf5/1.8.21
set_mod $1 pkg $siesta/Docs/build/netcdf/4.7.4
set_mod $1 pkg $siesta/Docs/build/zlib/1.2.11
set_mod $1 bin $siesta/Obj
unset siesta
reset_env $1 OMP_NUM_THREADS 1
```

加载设置时输入`bmod add siesta-v4.1-b4`，卸载时输入`bmod rm siesta-v4.1-b4`即可。bmod会自动把`$1`替换成`add`或`rm`。需注意，如果系统中安装了某个库的多个版本，而待运行的程序只能使用特定版本，那么在加载该版本前必须把其它版本卸掉。假设系统中还安装了hdf5-1.10.12，上述设置就要改成：

```bash
set_mod rm pkg $HOME/soft/lib/hdf5-1.10.12
siesta=$HOME/soft/dft/siesta-v4.1-b4
set_mod $1 pkg $siesta/Docs/build/flook/0.8.1
set_mod $1 pkg $siesta/Docs/build/hdf5/1.8.21
set_mod $1 pkg $siesta/Docs/build/netcdf/4.7.4
set_mod $1 pkg $siesta/Docs/build/zlib/1.2.11
set_mod $1 bin $siesta/Obj
unset siesta
reset_env $1 OMP_NUM_THREADS 1
```

由于我们希望在加载hdf5-1.8.21之前卸载hdf5-1.10.12，所以`set_mod`后面的参数必须是`rm`，不能用`$1`代替。

除了`add/rm`外，bmod命令支持如下选项：

- av: 列出所有可用脚本
- ls: 列出所有已加载脚本
- cl: 卸载已加载的脚本
- pg: 卸载所有脚本（含未加载）

`examples`目录下面中自带了很多用作例子的脚本。为了演示bmod的功能，我们把该目录临时加入`BMOD_MOD`环境变量：

```bash
set_env add BMOD_MOD $HOME/soft/bmod/examples 
```

然后运行`bmod av`命令，看一下输出结果：

```bash
[yhli@linux-h149 ~]$ bmod av
---- /home/yhli/soft/bmod/examples ----
   1) bgw-2.1
   2) boost-1.75.0
   3) fleurMaXR3.1
   4) office
   5) openmpi-4.0.5
   6) qe-6.6
   7) siesta-v4.1-b4
   8) spex05.00
   9) vesta
  10) vtst
  11) wannier90-2.1.0
  12) xcrysden-1.6.2
---- /home/yhli/soft/bmod/modules ----
   1) base
```

在输出中，给出了脚本所在的目录和该目录下的所有脚本。若有脚本保存在别的地方，可以用同样方法让bmod识别。

`bmod ls`命令的输出结果：

```bash
[yhli@linux-h149 ~]$ bmod ls
Currently loaded modules:
```

由于在一开始我们没有加载任何脚本，所以`bmod ls`输出为空。我们先用`bmod add`命令加载几个脚本后，再来看有何不同：

```bash
[yhli@linux-h149 ~]$ bmod add bgw-2.1 boost vesta
[yhli@linux-h149 ~]$ bmod ls
Currently loaded modules:
   1) bgw-2.1
   2) boost-1.75.0
   3) vesta
```

可以看到这几个脚本已被加载了。如果脚本名字中的版本号和主体以-隔开，可以省略不写。所以`bmod add boost`是可以的，而`bmod add vasp`却不行，必须输入全称`bmod add vasp.5.4.4`。

我们用`bmod cl`命令卸载已加载的脚本，看下效果：

```bash
[yhli@linux-h149 ~]$ bmod ls
Currently loaded modules:
```

可见所有已加载的脚本都已被卸载。`bmod pg`有类似效果，就不再展示了。

## 4. 补充说明

### 4.1 常见环境变量及报错信息

#### 4.1.1 PATH

系统需要通过这个环境变量查找可执行程序（二进制文件或脚本）。若未正确设置，会报类似如下错误：

```bash
[yhli@linux-h149 ~]$ pw.x
-bash: pw.x: command not found
```

#### 4.1.2 LD\_RUN\_PATH, LIBRARY\_PATH

这两个环境变量指定编译程序时，所用到的库文件所在位置。库文件有两种：静态链接库(`libxxx.a`)和动态链接库(`libxxx.so`)。若未正确设置，会报类似如下错误：

```bash
/usr/lib64/gcc/x86_64-suse-linux/7/../../../../x86_64-suse-linux/bin/ld: cannot find -lboost_python38
collect2: error: ld returned 1 exit status
make: *** [Makefile:65：pyxaid_core.so] Error 1
```

#### 4.1.3 LD\_LIBRARY\_PATH

这个环境变量指定运行程序时，所用到的动态链接库所在位置。若未正确设置，会报类似如下错误：

```bash
[yhli@linux-h149 doublecmd]$ ./doublecmd
./doublecmd: error while loading shared libraries: libQt5Pas.so.1: cannot open shared object file: No such file or directory
```

#### 4.1.4 C\_INCLUDE\_PATH, CPLUS\_INCLUDE\_PATH

这两个环境变量指定编译程序时，所用到的头文件所在位置。若未正确设置，会报类似如下错误：

```bash
pyxaid_core.cpp:10:10: fatal error: boost/python.hpp: No such file or directory
 #include <boost/python.hpp>
          ^~~~~~~~~~~~~~~~~~
compilation terminated.
make: *** [Makefile:59：pyxaid_core.o] Error 1
```

### 4.2 注意事项

目前bmod缺失对已加载模块检查的功能。以下述命令为例：

```bash
[yhli@linux-h149 ~]$ bmod add qe
[yhli@linux-h149 ~]$ bmod add fleurMaXR3.1
[yhli@linux-h149 ~]$ bmod ls
Currently loaded modules:
  1) qe-6.6
  2) fleurMaXR3.1
```

当执行第一条命令`bmod add qe`时，bmod加载了与qe相关的设置。一般来说，这需要卸载冲突库和加载依赖库。当执行第二条命令`bmod add fleurMaXR3.1`时，bmod同样会执行类似操作。但是，如果fleur依赖的某个库恰巧与qe冲突，或者与fleur有冲突的某个库恰巧是qe的依赖库，那么执行完第二句后，qe的环境变量就被破坏了。也就是说，bmod目前只能“狗熊掰玉米”，保证最后一个加载的程序可用。要解决这个问题，需要复杂的依赖和冲突关系分析。这已经超出了bmod的设计初衷，如有这方面的需求，可以选择功能更强大的Environment Modules系统，或本人开发的Pmod。
