workspace "naiveFFT"
	architecture "x64"

	configurations
	{
		"Debug",
		"Release"
	}

	
outputdir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"


project "naiveFFT"

	kind "WindowedApp"
	language "C++"
	cppdialect "C++17"

	targetdir ("bin/" .. outputdir .. "/%{prj.name}")
	objdir("bin-int/" .. outputdir .. "/%{prj.name}")
	debugdir ("bin/" .. outputdir .. "/%{prj.name}")


	warnings "High"

	files
	{
		"%{prj.name}/src/**.h",
		"%{prj.name}/src/**.cpp",
		"%{prj.name}/include/**.h",
		"%{prj.name}/include/**.cpp"

	}

	includedirs
	{
		"%{prj.name}/include",
	}
		filter "system:windows"

		systemversion "latest"
		staticruntime "on"
		flags {"MultiProcessorCompile"}



		filter "configurations:Debug"
			symbols "on"

		filter "configurations:Release"
			optimize "on"
			