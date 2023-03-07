const Mat33{T} = SMatrix{3,3,T,9}
const Mat22{T} = SMatrix{2,2,T,4}

const Vec2{T} = SVector{2,T}
const Vec3{T} = SVector{3,T}
const Vec4{T} = SVector{4,T}
const Vec6{T} = SVector{6,T}

const ID3 = Diagonal(Vec3(1,1,1))
const Mat224{T} = SMatrix{2,24,T,48}
const Mat324{T} = SMatrix{3,24,T,72}


const Mat44{T} = SMatrix{4,4,T,16} 
const Mat66{T} = SMatrix{6,6,T,36} 
const Mat77{T} = SMatrix{7,7,T,49}
const Mat36{T} = SMatrix{3,6,T,18}
const Mat312{T} = SMatrix{3,12,T,36}
const Mat612{T} = SMatrix{6,12,T,72}
const Mat712{T} = SMatrix{7,12,T,84}
const Mat1212{T} = SMatrix{12,12,T,144}
const Mat123{T} = SMatrix{12,3,T,36}


const Vec7{T} = SVector{7,T}
const Vec12{T} = SVector{12,T} 
const Vec24{T} = SVector{24,T} 

const ID3 = Diagonal(Vec3(1,1,1))