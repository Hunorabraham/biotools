const std = @import("std");

const Sequence = struct{
    name : []const u8,
    letters : []const u8,
    current: usize,
    pub fn next(Self: *Sequence) ?u8{
        if(Self.current >= Self.letters.len) return null;
        if(Self.letters[Self.current] == '\n') Self.current += 1; //skip newlines
        if(Self.letters[Self.current] == '\r') Self.current += 1; //skip carriage returns
        if(Self.letters[Self.current] == 'U') return 'T'; //return timin instead of uracil 
        if(Self.letters[Self.current] == 'u') return 't'; 
        if(Self.letters[Self.current] == 'N') return '-'; //return '-' instead of other "undefined" characters
        if(Self.letters[Self.current] == 'n') return '-'; 
        if(Self.letters[Self.current] == 'X') return '-'; 
        if(Self.letters[Self.current] == 'x') return '-'; 
        defer Self.current += 1;
        return Self.letters[Self.current];
    }
    pub fn jump(Self: *Sequence, n: usize) void{
        var count: usize = 0;
        while(count < n){
            count += 1;
            Self.current += 1;
            if(Self.letters[Self.current] == '\n') Self.current += 1;
            if(Self.letters[Self.current] == '\r') Self.current += 1;
        }
    }
    pub fn nextBatch(Self: *Sequence, comptime n:usize, target: *[n]u8) void{
        for(target, 0..) |t, i|{
            _ = t;
            target[i] = Self.next() orelse '-';
        }
    } 
    pub fn start(Self: *Sequence) void{
        Self.current = 0;
    }
    pub fn shortName(Self: *const Sequence) []const u8{
        for(Self.name,0..) |l, i|{
            if(l == '[') return Self.name[0..i];
        }
        return Self.name;
    }
};

const CompData = struct{
    name: []const u8,
    synonym: u32,
    nonSynonym: u32,
    pub fn ratio(Self: *const CompData) f32{
        if(Self.synonym == 0) return 0.0;
        if(Self.nonSynonym == 0) return 0.0;
        return @as(f32, @floatFromInt(Self.nonSynonym)) / @as(f32, @floatFromInt(Self.synonym));
    }
};


pub fn main() !void{
    //allocator
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    
    //arg reading
    var args = try std.process.argsWithAllocator(allocator);
    _ = args.skip();
    var arg1: ?[]const u8 = null;
    var arg2: ?[]const u8 = null;
    var arg3: ?[]const u8 = null;
    while(args.next()) |arg|{
        if(compstr(arg, "-n")){
            arg1 = args.next();
        }
        else if(compstr(arg, "-p")){
            arg2 = args.next();
        }
        else if(compstr(arg, "-h")){
            std.debug.print(\\switch syntax: -<switch> <parameter>
            \\-h
            \\  no parameter, displays this message
            \\-n
            \\  <file name>, specifies the file containing the sequences as nucleotides 
            \\-p
            \\  <file name>, specifies the file containing the sequences as amino acids 
            \\(note: sequences must be in the same order in both files, reference must be the first sequence)
            , .{});
            return;
        }
        else if(compstr(arg, "-o")){
            arg3 = args.next();
        }
        //std.debug.print("{s}\n",.{arg});
    }
    if(arg1 == null){
        std.debug.print("Specify file of nucleotide sequences with -n", .{});
        return;
    }
    if(arg2 ==  null){
        std.debug.print("Specify file of amino acid sequences with -p", .{});
        return;
    }
    const filenameN = arg1.?;
    const filenameP = arg2.?;
    
    //file reading
    const nucleotide = try readFile(filenameN, allocator);
    const protein = try readFile(filenameP, allocator);
    
    //tokenisation
    var sequences = [_]Sequence{Sequence{.name = "0", .letters = "0", .current = 0}}**100;
    var translations = [_]Sequence{Sequence{.name = "0", .letters = "0", .current = 0}}**100;
    const count = fillSequence(&sequences, nucleotide);
    std.debug.print("{d} sequences found\n", .{count});
    _= fillSequence(&translations, protein);
    
    //for(translations[0..count], sequences[0..count]) |t, s|{
    //    std.debug.print("name: {s}\nsequence len: {d}, translation len: {d}\n", .{s.shortName(), s.letters.len, t.letters.len});
    //}
    
    //calculate
    var results = [_]CompData{CompData{.name = "N", .synonym = 0, .nonSynonym = 0}}**100;
    var referenceN = sequences[0];
    var referenceP = translations[0];
    for(1..count) |i|{
        results[i-1].name = sequences[i].shortName();
        countSNS(&referenceN, &referenceP, &sequences[i], &translations[i], &results[i-1]);
    }
    for(results[0..(count-1)]) |r|{
        std.debug.print("{s}\n  Synonymous(dS): {d}\n  Nonsysnonymous(dN): {d}\n  ratio(dN/dS): {d}\n", .{r.name, r.synonym, r.nonSynonym, r.ratio()});
    }
}

pub fn countSNS(refN: *Sequence, refP: *Sequence, compN: *Sequence, compP: *Sequence, target: *CompData) void{
    var refTrip = [_]u8{0}**3;
    var compTrip = [_]u8{0}**3;
    refP.start();
    refN.start();
    while(refP.next()) |p|{
        var mut: u32 = 0;
        refN.nextBatch(3, &refTrip);
        compN.nextBatch(3, &compTrip);
        for(refTrip, compTrip) |r, c|{
            if(c == '-') continue;
            if(r == '-') continue;
            if(r == c) continue;
            mut += 1;
        }
        const pcomp = compP.next().?;
        if(p == '-') continue;
        if(pcomp == '-') continue;
        if(pcomp == p) {target.synonym += mut;}
        else {target.nonSynonym += mut;}
    }
}

fn compstr(a: []const u8, b: []const u8) bool{
    if(a.len != b.len) return false;
    for(a, b) |l, r|{
        if(l != r) return false;
    }
    return true;
}

fn fillSequence(target: []Sequence, source: []const u8) usize{
    var count: usize = 0;
    var start: usize = 0;
    var i: usize = 0;
    
    while(i < source.len) : (count += 1){
        if(count == 100){
            std.debug.print("!!---only first 100 sequences read---!!", .{});
            break;
        }
        //name
        i += 1; //skip >
        start = i;
        while(source[i] != '\n'){i+=1;}//skip to new line
        target[count].name = source[start..i];
        i += 1; //go over new line
        if(source[i] == '\r') i+=1; //ignore carriage return (fuck u microsoft)
        //sequence
        start = i;
        var extra_spaces: u8 = 0;
        while(i < source.len){
            if(source[i] == '\n') extra_spaces += 1;
            if(source[i] == '\r') extra_spaces += 1;
            if(source[i] == '>') break;
            i+=1;
        }
        target[count].letters = source[start..(i-extra_spaces)];
    }
    return count;
}

fn readFile(filename: []const u8, allocator: std.mem.Allocator) ![]u8{
    const cwd = std.fs.cwd();
    var file = try cwd.openFile(filename, .{});
    defer file.close();
    const stat = try file.stat();
    const size = stat.size;
    const content = try allocator.alloc(u8, size);
    _ = try file.read(content);
    return content;
}