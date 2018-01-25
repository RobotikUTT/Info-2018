#ifndef WHATCHER_H
#define WATCHER_H

#include "../include/Module.h" //This class need to know the class above
#include "../include/Interface.h"



class Watcher{
  private:
    Watcher();


  public:

    static void CreateSingleton();
    ~Watcher();
    void Update();
};
extern Watcher* g_watcher_instance_ptr;

#endif
