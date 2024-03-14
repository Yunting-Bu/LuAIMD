# LuAIMD: An ab initio molecular dynamics program

LuAIMD的名字玩了个小双关，LuA是鲁A，也就是济南，正好连着AIMD...很冷

## 依赖与编译
`libcint`：https://github.com/sunqm/libcint

### 编译
1. 打开`make.sh`，将`/PATH/libcint.so`中的`PATH`改为实际路径。
2. `chmod +x make.sh`赋予权限。
3. 执行`./make.sh`进行编译。
4. 环境变量中，增加`export LUAIMD_HOME=/PATH/LuAIMD`，其中`PATH`为实际路径。
5. 为了能在任意路径通过指令`LuAIMD`执行程序，可以在环境变量中添加
    ```
    alias LuAIMD='/PATH/LuAIMD/LuAIMD'
    export PATH=$PATH:/PATH/LuAIMD
    ```
## 功能
### 量子化学
1. 支持HF与MP2(fc)，MP2(Full)，在输入文件中分别写为`hf  mp2(fc)  mp2(full)`。
2. 支持的基组见`basis`文件，可以直接去基组定义库拷贝，为Gaussian形式。
3. HF支持解析梯度，其余为数值梯度，写法为`analy  num`。
4. 支持Damp与DIIS加速。
5. 支持三种初猜：核哈密顿初猜`core`，GWH`GWH`以及拓展Huckel`Huckel`。
### 分子动力学
1. 支持速度Verlet算法`Velocity_Verlet`与蛙跳法`leapfrog`。
2. 支持Berendsen热浴`Berendsen`。
3. 输入文件均为fs。
### 其他
1. 支持计算指定的键长键角，输入`-1`不进行计算。

## 输入文件与输出文件
### 输入文件
以`H2O.inp`为例：
```
$QC_ctrl
    method:                  hf
    basis:                   6-31G**
    grad:                    analy
    guess:                   Huckel
    damp:                    0.5
    DIIS:                    .true.
$end 

$MD_ctrl
    method:                  Velocity_Verlet
    init_temp:               298.15
    dt:                      0.5
    Nstep:                   1000
    thermostat:              Berendsen
    bath_temp:               298.15
    con_time:                30.0
$end 

$MOLE_analy
    number_of_bond_length:   2
    list:                    1,2 1,3
    number_of_bond_angle:    1
    list:                    2,1,3
$end 

$geom
    0 1
    H2O.xyz
$end
```
需要提供`mol.xyz`文件并在同一文件夹下。
### 输出文件
会输出`mol.out`, `pos.xyz`, `mole_geom.dat`

`mol.out`: 输出电子结构与动力学信息。

`pos.xyz`: 输出运动过程的坐标。

`mole_geom.dat`: 输出键长键角信息，当`number_of_bond_length`与`number_of_bond_angle`都为`-1`时不产生此文件。
