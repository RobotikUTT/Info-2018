//#include "ros.h"
#include <../include/Interface.h>
#include <../include/Verification.h>
#include <../include/Module.h>
#include <../include/Watcher.h>
#include <unistd.h>
#include <time.h>
#include <iostream>

 //ros::NodeHandle nh;
int main(){

  Interface::CreateSingleton();
  Verification::CreateSingleton();
  Module::CreateSingleton();
  Watcher::CreateSingleton();

  while (1) {
    std::cout << "time: " << time(0) << std::endl;
    interface_instance_ptr->Update();
    verification_instance_ptr->Update();
    module_instance_ptr->Update();
    watcher_instance_ptr->Update();
    usleep(1000000);
  }
}
