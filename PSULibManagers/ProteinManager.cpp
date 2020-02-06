#include "ProteinManager.h"

ProteinManager* ProteinManager::instance = 0;

ProteinManager::ProteinManager()
{

}

ProteinManager* ProteinManager::getInstance()
{
	if (!instance)
		instance = new ProteinManager();
	return instance;
}

void ProteinManager::createConfigData(string path)
{
}
