/** ACA::Arena
haiyinpiao@qq.com
*/

#ifdef _MSC_VER
// We'll also define this to stop MSVC complaining about sprintf().
#define _CRT_SECURE_NO_WARNINGS
#pragma comment(lib, "Irrlicht.lib")
#endif

#include <irrlicht.h>
#include "driverChoice.h"
#include <cstdio>
#include <ctime>
#include <iostream>
#include "actor_grt_rtw/actor.h"

using namespace std;
using namespace irr;

// integration steps, in ms.
const float g_intstep=30;

float g_elpstime=0;

int MAXSTEPRELOAD=50;
int g_reloadstep=0;

const int ACTOR_NUM = 50;

class MyEventReceiver : public IEventReceiver
{
public:
	// This is the one method that we have to implement
	virtual bool OnEvent(const SEvent& event)
	{
		// Remember whether each key is down or up
		if (event.EventType == irr::EET_KEY_INPUT_EVENT)
			KeyIsDown[event.KeyInput.Key] = event.KeyInput.PressedDown;

		return false;
	}

	// This is used to check whether a key is being held down
	virtual bool IsKeyDown(EKEY_CODE keyCode) const
	{
		return KeyIsDown[keyCode];
	}
	
	MyEventReceiver()
	{
		for (u32 i=0; i<KEY_KEY_CODES_COUNT; ++i)
			KeyIsDown[i] = false;
	}

private:
	// We use this array to store the current state of each key
	bool KeyIsDown[KEY_KEY_CODES_COUNT];
};

int main(int argc, char *argv[])
{
    actorModelClass f16[ACTOR_NUM];

    for(int i=0; i<ACTOR_NUM; i++)
    {
        f16[i].initialize();

        // input
        f16[i].actor_U.nnk_c = float(rand()%5000)/5000.0f;
        f16[i].actor_U.nxk_c = float(rand()%1000)/1000.0f;
        f16[i].actor_U.mudot_c = float(rand()%700)/1000.0f;

        // initial condition
        f16[i].actor_P.xg0_Value[0] = rand()%7000;
        f16[i].actor_P.xg0_Value[1] = -20000+rand()%27000;
        f16[i].actor_P.xg0_Value[2] = -1*(rand()%10000);
        f16[i].actor_P.att_g0_Value[2] = (-180.0f+float(rand()%360))/57.3f;
    }


    // holds the simulation steps
    unsigned int k=0;

	MyEventReceiver receiver;

    IrrlichtDevice* device = createDevice(video::EDT_OPENGL,
            core::dimension2d<u32>(1024, 768), 16, false, false, false, &receiver);

	if (device == 0)
		return 1; // could not create selected driver.

	video::IVideoDriver* driver = device->getVideoDriver();
	scene::ISceneManager* smgr = device->getSceneManager();

    //data display fonts
    gui::IGUIFont* font = device->getGUIEnvironment()->getBuiltInFont();

    // add light
    scene::ISceneNode* plight_node = smgr->addLightSceneNode(0, core::vector3df(10000,10000,10000),
                                                             video::SColorf(1.0f, 0.6f, 0.7f, 1.0f), 1000000.0f);

    scene::IAnimatedMeshSceneNode* anms[ACTOR_NUM];

    for(int i=0; i<ACTOR_NUM; i++)
    {
        anms[i] =
            smgr->addAnimatedMeshSceneNode(smgr->getMesh("../aca_arena/T50-airframe.3ds"));

        if (anms[i])
        {
            anms[i]->setMaterialTexture(0, driver->getTexture("../aca_arena/T50-airframe.png"));
            /*scene::ISceneNodeAnimator* anim =
                smgr->createFlyCircleAnimator(core::vector3df(0,0,1), 10.0f);
            if (anim)
            {
                anms->addAnimator(anim);
                anim->drop();
            }*/

            anms[i]->setMaterialFlag(video::EMF_LIGHTING, false);

            anms[i]->setFrameLoop(0, 13);
            anms[i]->setAnimationSpeed(15);

            anms[i]->setScale(core::vector3df(2000.f,2000.f,2000.f));
            /*anms->setRotation(core::vector3df(0,0,0));
            anms->setPosition(core::vector3df(0,0,0));*/
    //		anms->setMaterialTexture(0, driver->getTexture("../../media/sydney.bmp"));

        }
    }

    scene::ITerrainSceneNode* terrain = smgr->addTerrainSceneNode(
                "../irrlicht-1.8.4/media/terrain-heightmap.bmp",
                0,                  // parent node
                -1,                 // node id
                core::vector3df(-40000.f, 0.f, -20000.f),     // position
                core::vector3df(90.f, 0.f, 0.f),     // rotation
                core::vector3df(300.f, 20.4f, 300.f),  // scale
                video::SColor ( 180, 180, 180, 180 ),   // vertexColor
                5,                  // maxLOD
                scene::ETPS_17,             // patchSize
                4                   // smoothFactor
                );

    terrain->setMaterialFlag(video::EMF_LIGHTING, false);
    /*terrain->setMaterialTexture(0,
                                driver->getTexture("../irrlicht-1.8.4/media/terrain-texture.jpg"));*/
    terrain->setMaterialTexture(1,
                                driver->getTexture("../irrlicht-1.8.4/media/detailmap3.jpg"));
    terrain->setMaterialType(video::EMT_DETAIL_MAP);
    terrain->scaleTexture(1.0f, 20.0f);

    terrain->setMaterialFlag(video::EMF_WIREFRAME,
            !terrain->getMaterial(0).Wireframe);
    terrain->setMaterialFlag(video::EMF_POINTCLOUD, false);

	/*
	To be able to look at and move around in this scene, we create a first
	person shooter style camera and make the mouse cursor invisible.
	*/
    //scene::ICameraSceneNode* camera=smgr->addCameraSceneNodeFPS();
    scene::ICameraSceneNode* camera=smgr->addCameraSceneNode(0, core::vector3df(30000,-0,-6000), core::vector3df(0,0,0));
    camera->setUpVector( core::vector3df(0,0,-1) );
    //camera->setRotation(core::vector3df(0,0,0.0f));
    device->getCursorControl()->setVisible(false);
    camera->setFarValue(142000.0f);

	int lastFPS = -1;

	// In order to do framerate independent movement, we have to know
	// how long it was since the last frame
	u32 then = device->getTimer()->getTime();

	// This is the movemen speed in units per second.
    const f32 MOVEMENT_SPEED = 100.f;

	while(device->run())
	{
		// Work out a frame delta time.
		const u32 now = device->getTimer()->getTime();
		const f32 frameDeltaTime = (f32)(now - then) / 1000.f; // Time in seconds
		then = now;

        g_elpstime += frameDeltaTime;


        if( g_elpstime>(0.003f/*g_intstep/10000.0f*/) )
        {
            for(int i=0; i<ACTOR_NUM; i++)
            {
                //output
                float x=f16[i].actor_Y.Xg[0];
                float y=f16[i].actor_Y.Xg[1];
                float z=f16[i].actor_Y.Xg[2];

                float mu=f16[i].actor_Y.att_g[0];
                float gamma=f16[i].actor_Y.att_g[1];
                float psi=f16[i].actor_Y.att_g[2];

                f16[i].step();

                anms[i]->setRotation(core::vector3df(mu*57.3+180,gamma*57.3,psi*57.3));
                anms[i]->setPosition(core::vector3df(x+5000,y+5000,z-2000));
            }

            g_elpstime=0;
            g_reloadstep++;
        }

        //re-initialize
        if( g_reloadstep>MAXSTEPRELOAD )
        {
            for(int i=0; i<ACTOR_NUM; i++)
            {
                f16[i].initialize();

                // input
                f16[i].actor_U.nnk_c = float(rand()%5000)/5000.0f;
                f16[i].actor_U.nxk_c = float(rand()%1000)/1000.0f;
                f16[i].actor_U.mudot_c = float(rand()%700)/1000.0f;

                // initial condition
                f16[i].actor_P.xg0_Value[0] = rand()%7000;
                f16[i].actor_P.xg0_Value[1] = -20000+rand()%27000;
                f16[i].actor_P.xg0_Value[2] = -1*(rand()%10000);
                f16[i].actor_P.att_g0_Value[2] = (-180.0f+float(rand()%360))/57.3f;
            }

            g_reloadstep=0;
        }

        driver->beginScene(true, true, video::SColor(80,80,80,80));

		smgr->drawAll(); // draw the 3d scene
		device->getGUIEnvironment()->drawAll(); // draw the gui environment (the logo)

        // draw flight parameters text.
        if (font)
        {
            font->draw(L"This demo shows a hyper real-time f-16 dynamics simulation.\n haiyinpiao@qq.com.\n %f %f %f, x,y,z\n %f %f %f, x,y,z\n %f %f %f, x,y,z\n %f %f %f, x,y,z",
                core::rect<s32>(100,100,800,400),
                video::SColor(255,255,255,255));
        }

		driver->endScene();

		int fps = driver->getFPS();

		if (lastFPS != fps)
		{
            core::stringw tmp(L"ACA::Arena - haiyinpiao@qq.com [");
			tmp += driver->getName();
			tmp += L"] fps: ";
			tmp += fps;

			device->setWindowCaption(tmp.c_str());
			lastFPS = fps;
		}
	}

	/*
	In the end, delete the Irrlicht device.
	*/
	device->drop();

    for(int i=0; i<ACTOR_NUM; i++)
    {
        f16[i].terminate();
    }
	
	return 0;
}

/*
That's it. Compile and play around with the program.
**/
