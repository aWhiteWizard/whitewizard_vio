# whitewizard_vio
a simple vio system for practice

//  运行效果展示

  链接: https://pan.baidu.com/s/1bkSIZfHv7A_nhNgOSV7M3w?pwd=slam 提取码: slam 

//  下载与编译

  git clone https://github.com/aWhiteWizard/whitewizard_vio.git 
  
  mkdir build && cd build
  
  cmake ..
  
  make
  
//  运行方法

  这套代码作为个人练习成果，没有做优化，只是实现了功能，崩溃很多，一次起不起来多启动几次，总会有一次成功的
  
  cd whitewizard_vio/build
  
  ./test_euroc euroc数据集路径 ../config/
 
  example:
  
  ./test_euroc /home/whitewizard/4_dataset/MH_05_difficult/mav0/ ../config/


//  结语

  终于，历时两个月，基于vins-mono从头到尾手动实现了一套vio系统，虽然LM阻尼因子收敛仍然有问题，偶尔还是会从几万收敛到几万，虽然偶尔会初始化不过，虽然长时间不更新地图点会偶尔崩溃，虽然可能是由于划窗用舒尔补可能是数据有积累误差导致频繁重启，sfm偶尔也会实效。。。。

  即便这么多的bug还没处理，但是总之，它还是踉踉跄跄的跑起来了，算是圆了自己把vio从头到尾的实现一遍的梦想。也要感谢在高博和贺博的帮助，以及感慨以下vins-mono代码写的确实不错，如果真的让我自己从头到尾弄帕是一年也搞不完。

  之后对vio的研究还是基于vins-mono吧，之前研究了一下rootba并且写了demo对照ceres和rootba的速度，发现rootba提升很大，可以研究一下提升vins-mono的效率

//基础源代码连接
  
  VINS-Mono https://github.com/HKUST-Aerial-Robotics/VINS-Mono
 
