#include <winsock2.h>

int Startup()
{
    WSADATA wsadata;
	return (WSAStartup(MAKEWORD(1,1), &wsadata) == SOCKET_ERROR);
}

int Cleanup()
{
    return WSACleanup();
}

char *getUserName(char *user, int len)
{
    DWORD bsize = len;
	return GetUserName(user,&bsize) ? user : "Windows User";
}
