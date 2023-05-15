# Typora Linux破解

## 1.下载自动注入工具与破解验证工具

+ 自动注入工具：[自动注入工具][https://github.com/DiamondHunters/NodeInject]

+ 破解验证脚本：[破解验证脚本][https://github.com/DiamondHunters/NodeInject_Hook_example]

## 2.替换脚本

将破解验证脚本中的`hook.js`改名`hooklog.js`并替换到自动注入工具的src中

## 3.安装cargo

```sh
curl https://sh.rustup.rs -sSf | sh
```

## 4.编译自动注入脚本

在`NodeInject`目录下运行`cargo build --features no_embed`构建工具

## 5. 运行工具

将`NodeInject`目录下所有文件拷贝到typora根目录`/usr/share/typora`并运行自动注入工具`cargo run`，运行成功后在`NodeInject_Hook_example/license-gen`目录下运行`cargo run`生成激活码，用此激活码激活typora